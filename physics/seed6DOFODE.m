function [xdot, intermediates] = seed6DOFODE(t, x, seedParams)
% SEED6DOFODE  Right-hand side of the 3D seed dynamics ODE (for ode45).
%
% Strip-theory quasi-steady aerodynamics extending the 2D Andersen-Pesavento-Wang
% model to 6-DOF (13-state) flight. Computes the state derivative; optionally
% returns per-strip intermediate quantities for plotting/debugging (ode45 ignores
% the second output during integration, so re-call this on the saved states in
% post-processing to recover them).
%
% STATE  x = [ r(3) ; q(4) ; v(3) ; omega(3) ]  (13x1)
%   r     : CoM position, inertial frame (m)
%   q     : orientation quaternion [q0;q1;q2;q3], body->world, scalar-first
%   v     : CoM velocity, inertial frame (m/s)
%   omega : angular velocity, BODY frame (rad/s)
% NOTE: r/v are the CoM in the inertial frame; mp.c (below) is the CoM location
% WITHIN the body (body-datum coords) used for moment arms -- a different quantity.
%
% INPUTS
%   t          : time (s).
%   x          : 13x1 state vector (above).
%   seedParams : struct from setupSeedShapeAndMass, plus the fields:
%                  .rhoFluid  fluid density (kg/m^3), e.g. 1.225 for air
%                  .g         gravity magnitude (m/s^2), e.g. 9.81
%                  .aero      (optional) empirical coefficient struct for
%                             computeAeroCoeffs; defaults used if absent.
%
% OUTPUTS
%   xdot          : 13x1 state derivative.
%   intermediates : struct of per-strip and total quantities (see end of file).

% =========================================================================
% 0. PARSE STATE
% =========================================================================
r     = x(1:3);        %#ok<NASGU>  % CoM position (not needed for the derivative)
q     = x(4:7);
v     = x(8:10);       % CoM velocity, inertial
omega = x(11:13);      % angular velocity, body frame

q = q / norm(q);       % normalise for the physics (state may have drifted)
R = quatToRotm(q);     % body->world rotation

% Aerodynamic coefficient parameters (default inside computeAeroCoeffs if absent).
if isfield(seedParams, 'aero'); aeroParams = seedParams.aero; else; aeroParams = []; end
rhoFluid = seedParams.rhoFluid;
gMag     = seedParams.g;

% =========================================================================
% 1. MASS + ADDED MASS at this time
% =========================================================================
mp = getMassProperties(t, seedParams);   % mp.M, mp.c, mp.I_G, mp.I_G_dot, mp.A_trans, mp.A_rot

% =========================================================================
% 2. KINEMATICS: body-frame velocity at each strip
% =========================================================================
% Transports the CoM velocity to each strip via the transport theorem; returns
% body-frame velocities (3xM). We use the total velocity (its in-plane part) for
% the aerodynamics.
localVels = computeSeedLocalVel(seedParams, mp.c, q, v, omega);
stripVel  = localVels.totalVel;           % 3xM, body frame

% =========================================================================
% 3. STRIP LOOP: forces and torques about the CoM (body frame)
% =========================================================================
chord = seedParams.strips.chord;
dz    = seedParams.strips.dz;
xGeo  = seedParams.strips.xgc_body;       % strip geometric-centre chordwise position (body x)
zGeo  = seedParams.strips.zgc_body;       % strip geometric-centre spanwise position (body z)
liftMult = seedParams.strips.liftMult;    % 1xM per-strip translational-lift multipliers
dragMult = seedParams.strips.dragMult;    % 1xM per-strip drag multipliers
numStrips = numel(chord);

F_aero_body = [0; 0; 0];
tau_body    = [0; 0; 0];

% Pre-allocate intermediate storage.
alphaAll  = zeros(1, numStrips);
xcpAll    = zeros(1, numStrips);
CT_all    = zeros(1, numStrips);
CD_all    = zeros(1, numStrips);
Fbody_all = zeros(3, numStrips);
Tbody_all = zeros(3, numStrips);

omega_z = omega(3);   % only spanwise spin drives the 2D rotational lift/damping

% =========================================================================
% 3a. WHOLE-SEED SPIN-DAMPING ABOUT THE CHORDWISE AND NORMAL AXES
%     (x-axis and y-axis analogues of Tr; Tr itself only damps omega_z)
% =========================================================================
% stripSpinDamping (2D Tr) only damps spin about the SPANWISE axis (omega_z).
% Rotation about the CHORDWISE axis (omega_x) induces a velocity that varies
% with SPANWISE position instead, so the analogous quartic integral is redone
% over the span (spanSpinDamping.m). Rotation about the NORMAL axis (omega_y)
% has no clean analogous closed form (see normalSpinDamping.m) and uses a
% placeholder quadratic form instead. Both are evaluated ONCE for the whole
% seed rather than per strip -- see spanSpinDamping.m for why.
omega_x = omega(1);   % chordwise angular velocity (rad/s)
omega_y = omega(2);   % plate-normal angular velocity (rad/s)

totalSpan = sum(dz);               % total spanwise extent of the seed (m)
wingArea  = sum(chord .* dz);      % total wing area (m^2)
chordMean = wingArea / totalSpan;  % mean aerodynamic chord (m); exact for a
                                    % non-tapered wing, representative otherwise
zGeoCenterSpan = sum(chord .* dz .* zGeo) / wingArea;  % area-weighted spanwise
                                                        % geometric centre (m)
xGeoCenterChord = sum(chord .* dz .* xGeo) / wingArea; % area-weighted chordwise
                                                        % geometric centre (m)
comSpanOffset  = mp.c(3) - zGeoCenterSpan;   % CoM spanwise offset from that centre (m)

% CD_rot, CD0, and C_fy are all alpha-independent (see computeAeroCoeffs), so
% any alpha works here; fetch all three from a single call.
constCoeffs  = computeAeroCoeffs(0, aeroParams);
CD_rot_const = constCoeffs.CD_rot;
CD0_const    = constCoeffs.CD0;
C_fy_const   = constCoeffs.C_fy;

Tx = spanSpinDamping(totalSpan, chordMean, omega_x, comSpanOffset, CD_rot_const, rhoFluid);

% Characteristic radius for the normal-axis placeholder: half the total span (for now).
R_normalSpin = totalSpan / 2;
Ty = normalSpinDamping(R_normalSpin, omega_y, CD0_const, C_fy_const, rhoFluid);

% =========================================================================
% 3b. WHOLE-SEED SPANWISE-FLOW FORCE  (optional; fills the unmodeled body-z force)
% =========================================================================
% The chordwise strips produce zero body-z (spanwise) force, so a sideways-
% sliding seed feels no aerodynamic resistance. When enabled, computeSpanForce
% supplies a single whole-seed force from flow in the span-normal (z-y) plane.
% Default OFF, so existing (validated, in-plane) behaviour is unchanged unless
% explicitly requested via seedParams.enableSpanForce.
enableSpanForce = isfield(seedParams, 'enableSpanForce') && seedParams.enableSpanForce;

% Contributions default to zero so the post-loop accumulation is unconditional.
F_span_apply = [0; 0; 0];   % force  contribution (added after the strip loop)
tau_span     = [0; 0; 0];   % torque contribution (added after the strip loop)
F_span_full  = [0; 0; 0];   % full span force, for intermediates (zero when disabled)

if enableSpanForce
    % Whole-seed geometric centre in body coords (area-weighted; body y = 0).
    seedGeoCenter = [xGeoCenterChord; 0; zGeoCenterSpan];

    % Bulk velocity transported to the geometric centre (transport theorem):
    %   v_gc^body = R' * v_com^inertial + omega^body x (geoCentre - CoM)
    v_com_body = R.' * v;
    v_gc_body  = v_com_body + cross(omega, seedGeoCenter - mp.c);
    vSpan_gc   = v_gc_body(3);   % spanwise    (body z) component
    vNormal_gc = v_gc_body(2);   % plate-normal (body y) component

    % Angle of attack in the span-normal plane: atan2(vNormal, vSpan). Reuse
    % computeAngleOfAttack by placing span in its "chord" (first) slot.
    beta       = computeAngleOfAttack([vSpan_gc; vNormal_gc; 0]);
    spanCoeffs = computeAeroCoeffs(beta, aeroParams);

    % Full span force [0; Fy; Fz] and its span-CoP fraction (fraction of span).
    [F_span_full, l_cp_frac_span] = computeSpanForce(vSpan_gc, vNormal_gc, ...
        totalSpan, chordMean, spanCoeffs.C_span, spanCoeffs, rhoFluid);

    % --- FORCE: keep ONLY the body-z (spanwise) component ------------------
    % Explicitly DISCARD body-x (already 0) and body-y (normal). The strips
    % already model the normal-direction force; re-adding F_span_full(2) would
    % double-count it. Only the body-z direction was previously unmodeled.
    F_span_apply = [0; 0; F_span_full(3)];

    % --- TORQUE: use the FULL force, including the discarded body-y --------
    % Deliberate Option-B choice: F_span_full(2) is NOT summed into the net
    % force above, but it IS used for the moment arm here. The strips supply the
    % normal force and its chordwise-CoP moment, but strip theory assumes
    % uniform spanwise loading and so can never produce the ROLL moment (about
    % body x) from a span-offset pressure centre; crossing the migrating span-
    % CoP arm with the full force recovers exactly that missing roll moment.
    % (Were we to use only F_span_apply's z-component, the span-CoP -- which
    % migrates along body z, parallel to that force -- would add no torque, and
    % this would collapse to a pure yaw from any chordwise CoM offset, i.e. the
    % same as applying at the geometric centre.)
    zSpanCoP  = computeStripCoP(l_cp_frac_span, totalSpan, zGeoCenterSpan);  % migrating span CoP (body z)
    r_spanCoP = [xGeoCenterChord; 0; zSpanCoP] - mp.c;   % span-CoP arm, relative to CoM
    % r_spanCoP = seedGeoCenter - mp.c;   % <-- ALTERNATIVE: apply the torque at the geometric centre instead
    tau_span  = cross(r_spanCoP, F_span_full);
end

for i = 1 : numStrips

    % --- strip velocity (body frame); discard spanwise component for the 2D aero
    vStrip  = stripVel(:, i);
    vChord  = vStrip(1);
    vNormal = vStrip(2);

    % --- angle of attack
    alpha = computeAngleOfAttack(vStrip);

    % --- aerodynamic coefficients at this angle of attack
    coeffs = computeAeroCoeffs(alpha, aeroParams);

    % --- dynamic centre of pressure (chordwise), and the geometric-centre reference
    xcp         = computeStripCoP(coeffs.l_cp_frac, chord(i), xGeo(i));
    r_cp        = [xcp;     0; zGeo(i)] - mp.c;   % CoP relative to CoM (arm for lift+drag)
    r_geoCenter = [xGeo(i); 0; zGeo(i)] - mp.c;   % geometric centre rel. to CoM (arm for rot. lift)

    % --- forces (body frame): translational lift+drag, and rotational lift
    % Per-strip liftMult/dragMult scale this strip's translational lift and drag
    % only (not the rotational lift).
    [F_transl, F_rotLift] = computeStripForces(vChord, vNormal, chord(i), dz(i), ...
                                               omega_z, coeffs, rhoFluid, ...
                                               liftMult(i), dragMult(i));
    dF = F_transl + F_rotLift;

    % --- torques about the CoM: each force at its own chordwise point
    tau_transl  = cross(r_cp,        F_transl);   % 2D Tt   (at CoP)
    tau_rotLift = cross(r_geoCenter, F_rotLift);  % 2D Tcr  (at geometric centre / mid-chord)

    % --- extra torque: nonlinear spin damping about the spanwise axis (2D Tr)
    comOffsetFromGeoCenter = mp.c(1) - xGeo(i);   % CoM chordwise offset from this strip's geometric centre
    Tr        = stripSpinDamping(chord(i), dz(i), omega_z, comOffsetFromGeoCenter, coeffs.CD_rot, rhoFluid);
    tau_spin  = [0; 0; Tr];

    dTau = tau_transl + tau_rotLift + tau_spin;

    % --- accumulate
    F_aero_body = F_aero_body + dF;
    tau_body    = tau_body    + dTau;

    % --- store intermediates
    alphaAll(i)   = alpha;
    xcpAll(i)     = xcp;
    CT_all(i)     = coeffs.CT;
    CD_all(i)     = coeffs.CD;
    Fbody_all(:,i)= dF;
    Tbody_all(:,i)= dTau;
end

% Add the whole-seed contributions (sections 3a, 3b) -- single contributions
% for the whole seed, not per-strip ones, so they are added here rather than
% inside the strip-accumulation loop above:
%   - spin-damping torques Tx (chordwise) and Ty (normal)
%   - the optional spanwise-flow force (body-z only) and its torque
tau_body    = tau_body    + [0; Ty; 0] + tau_span; %ZEROED tx
F_aero_body = F_aero_body + F_span_apply;

% =========================================================================
% 4. SUM FORCES (inertial frame) + GRAVITY  ->  linear acceleration
% =========================================================================
% Aero forces were summed in the body frame; rotate the total to inertial.
F_aero_inertial = R * F_aero_body;

% Gravity along world -Y (Y-up convention). Buoyancy neglected for seeds.
g_world = [0; -gMag; 0];
F_grav  = mp.M * g_world;

F_total_inertial = F_aero_inertial + F_grav;

a_inertial = translationDynamics(F_total_inertial, R, mp);

% =========================================================================
% 5. TORQUES (body frame)  ->  angular acceleration
% =========================================================================
% Torques stay in the body frame because omega and alpha are body-frame; the
% modified Euler equation handles the gyroscopic and I_G_dot terms.
alpha_body = rotationDynamics(tau_body, omega, mp);

% =========================================================================
% 6. QUATERNION KINEMATICS  ->  orientation derivative
% =========================================================================
dq = quatKinematics(q, omega);

% =========================================================================
% 7. ASSEMBLE STATE DERIVATIVE
% =========================================================================
% dr/dt = v (inertial velocity is a state, so this is direct).
xdot = [ v ; dq ; a_inertial ; alpha_body ];

% =========================================================================
% 8. OPTIONAL INTERMEDIATES (for post-processing / plotting)
% =========================================================================
if nargout > 1
    intermediates.alpha           = alphaAll;         % 1xM angle of attack (rad)
    intermediates.xcp             = xcpAll;           % 1xM chordwise CoP (m)
    intermediates.CT              = CT_all;           % 1xM lift coefficient
    intermediates.CD              = CD_all;           % 1xM drag coefficient
    intermediates.F_strip_body    = Fbody_all;        % 3xM strip force, body frame
    intermediates.tau_strip_body  = Tbody_all;        % 3xM strip torque, body frame
    intermediates.F_aero_body     = F_aero_body;      % 3x1 total aero force, body
    intermediates.F_aero_inertial = F_aero_inertial;  % 3x1 total aero force, inertial
    intermediates.F_total_inertial= F_total_inertial; % 3x1 aero + gravity, inertial
    intermediates.Tx_spanSpin     = Tx;               % whole-seed chordwise spin-damping torque (N*m)
    intermediates.Ty_normalSpin   = Ty;               % whole-seed normal-axis spin-damping torque (N*m)
    intermediates.F_span_full     = F_span_full;      % 3x1 full span-flow force, body (y-component discarded from sum)
    intermediates.F_span_apply    = F_span_apply;     % 3x1 span-flow force actually added (body-z only)
    intermediates.tau_span        = tau_span;         % 3x1 span-flow torque, body (uses full force + migrating span CoP)
    intermediates.tau_body        = tau_body;         % 3x1 total torque, body
    intermediates.a_inertial      = a_inertial;       % 3x1 linear acceleration
    intermediates.alpha_body      = alpha_body;       % 3x1 angular acceleration
end

end
