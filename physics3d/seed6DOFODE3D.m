function [xdot, intermediates] = seed6DOFODE3D(t, x, seedParams)
% SEED6DOFODE3D  Right-hand side of the NON-PLANAR (3D-shape) seed dynamics ODE.
%
% Extension of the planar seed6DOFODE to a seed whose strips are NOT confined to
% the body x-z plane: each strip carries a 3D geometric-centre position AND its
% own orthonormal local frame (chord/normal/span directions) in body coordinates,
% so twist / camber / dihedral can be represented. The strip aerodynamics are the
% SAME 2D coefficient laws as the planar model, evaluated in each strip's local
% (chord, normal) plane and then rotated back to the body frame; the net body
% force/torque go to the shared rigid-body core (rigidBody6DOF), exactly as in the
% planar RHS.
%
% SCAFFOLD STATUS: this reduces EXACTLY to the planar model when every strip's
% frame is identity and ygc = 0 (i.e. a flat plate) -- the flat-equivalence
% regression. Twist/curvature are introduced by the builder (setupSeedShape3D)
% populating non-identity frames; the whole-seed span-force block below is planar
% heritage (valid only for a flat seed) and should be switched OFF
% (enableSpanForce = false) for genuinely curved seeds, which produce their
% spanwise force from geometry instead.
%
% Requires the 'shape3d' seedParams contract (see validateSeedParams): the planar
% strip fields PLUS strips.ygc_body and strips.chordDir/normalDir/spanDir. If a
% planar seed (no frames) is passed, identity frames + y=0 are assumed, so this
% function also reproduces the planar model on a planar seed.
%
% Optionally returns per-strip intermediate quantities for plotting/debugging
% (ode45 ignores the second output during integration).
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

% --- 3D strip geometry: out-of-plane position + per-strip local frame -----
% ygc is the strip centre's body-y (0 for a flat plate). chordDir/normalDir/
% spanDir are the strip's local axes expressed in BODY coords (columns of the
% strip->body rotation). A flat plate has ygc=0 and the identity frame; a planar
% seedParams that lacks these fields is treated as exactly that, so this RHS
% reproduces the planar model on a planar seed.
if isfield(seedParams.strips, 'ygc_body')
    yGeo = seedParams.strips.ygc_body;
else
    yGeo = zeros(1, numStrips);
end
if isfield(seedParams.strips, 'chordDir')
    chordDir  = seedParams.strips.chordDir;    % 3xM
    normalDir = seedParams.strips.normalDir;   % 3xM
    spanDir   = seedParams.strips.spanDir;     % 3xM
else
    chordDir  = repmat([1;0;0], 1, numStrips);
    normalDir = repmat([0;1;0], 1, numStrips);
    spanDir   = repmat([0;0;1], 1, numStrips);
end

F_aero_body = [0; 0; 0];
tau_body    = [0; 0; 0];

% Pre-allocate intermediate storage.
alphaAll  = zeros(1, numStrips);
xcpAll    = zeros(1, numStrips);
CT_all    = zeros(1, numStrips);
CD_all    = zeros(1, numStrips);
Fbody_all = zeros(3, numStrips);
Tbody_all = zeros(3, numStrips);

% (The 2D rotational lift/damping is driven by spin about each strip's LOCAL span
% axis; that is computed per strip inside the loop as spanDir(:,i)'*omega. For a
% flat plate spanDir = [0;0;1], recovering the planar omega_z = omega(3).)

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
C_Tx_const   = constCoeffs.C_Tx;

Tx = spanSpinDamping(totalSpan, chordMean, omega_x, comSpanOffset, CD_rot_const, rhoFluid);

% Characteristic radius for the normal-axis placeholder: half the total span (for now).
R_normalSpin = totalSpan / 2;
Ty = normalSpinDamping(R_normalSpin, omega_y, CD0_const, C_fy_const, rhoFluid);

% =========================================================================
% 3b. WHOLE-SEED SPANWISE-FLOW FORCE  (optional; fills the unmodeled body-z force)
% =========================================================================
% The chordwise strips produce zero body-z (spanwise) force, so a sideways-
% sliding seed feels no aerodynamic resistance. computeSpanForce supplies a
% single whole-seed force from flow in the span-normal (z-y) plane.
%
% Three independent switches (all DEFAULT TRUE when the field is absent):
%   .enableSpanForce        - apply the span force at all
%   .enableSpanGeomVelocity - include the rotational omega x r transport when
%                             sampling the velocity at the geometric centre.
%     Set FALSE to drive the span force from the CoM translational velocity
%     only. Diagnostic value: the omega x r term feeds a SPANWISE velocity only
%     when the geometric-centre-to-CoM arm has a CHORDWISE (or out-of-plane)
%     component, since (omega x r)_z = omega_x*r_y - omega_y*r_x. For a purely
%     spanwise CoM offset it contributes nothing to this force, so toggling
%     this isolates how much the rotational sweep (vs. the descent velocity
%     resolved in the spinning body frame) is driving the span force.
%   .enableSpanCOPMigration - place the span force at the MIGRATING span centre
%     of pressure (true) or at the fixed geometric centre (false). The migrating
%     CoP arm scales with the TOTAL SPAN and its l_cp_frac(beta) oscillates at
%     the spin frequency during autorotation, which can parametrically
%     destabilise the roll axis; setting this false removes that oscillating
%     arm while leaving the span FORCE untouched.
%
% The span torque is additionally scaled by the aero constant C_span_torque
% (default 1), so the force and torque contributions can be tuned INDEPENDENTLY
% -- useful because tumbling needs the force while autorotation is sensitive to
% the torque.
% Physics switches -- absent fields fall back to the "FULL-minus-geomVelocity"
% default (span force + torque + CoP migration + Tx on; geom-velocity and
% attenuation off), so the model always defaults to that config.
if isfield(seedParams, 'enableSpanForce')
    enableSpanForce = seedParams.enableSpanForce;
else
    enableSpanForce = true;
end
if isfield(seedParams, 'enableSpanGeomVelocity')
    useSpanGeomVelocity = seedParams.enableSpanGeomVelocity;
else
    useSpanGeomVelocity = false;
end
if isfield(seedParams, 'enableSpanCOPMigration')
    useSpanCOPMigration = seedParams.enableSpanCOPMigration;
else
    useSpanCOPMigration = true;
end
% enableSpanTorque (default true when the span force is on): apply the span
% force's TORQUE contribution. Set FALSE to keep the span FORCE (the body-z
% resistance) but drop its moment entirely -- lets you test whether the force
% alone is enough and the torque is what destabilises autorotation.
if isfield(seedParams, 'enableSpanTorque')
    useSpanTorque = seedParams.enableSpanTorque;
else
    useSpanTorque = true;
end
% enableSpanTorqueAttenuation (default false): scale the span torque by a
% reduced-frequency factor 1/(1+(k/k0)^2), k = |omega_y|*S/(2*v_ip). Keys on the
% NORMAL-axis spin omega_y (the autorotation spin) rather than total |omega|, so
% it suppresses the span torque during fast autorotation -- where its spin-
% frequency oscillation parametrically destabilises the roll axis -- while
% leaving it intact during tumbling (an omega_z mode). Set FALSE for the raw
% (unattenuated) span torque.
if isfield(seedParams, 'enableSpanTorqueAttenuation')
    useSpanTorqueAttenuation = seedParams.enableSpanTorqueAttenuation;
else
    useSpanTorqueAttenuation = false;
end
% enableTxDamping (default TRUE): apply the Tx chordwise-axis (roll) spin-
% damping term (spanSpinDamping.m), scaled by aero.C_Tx. On by default.
if isfield(seedParams, 'enableTxDamping')
    useTxDamping = seedParams.enableTxDamping;
else
    useTxDamping = true;
end

% Contributions default to zero so the post-loop accumulation is unconditional.
F_span_apply    = [0; 0; 0];   % force  contribution (added after the strip loop)
tau_span        = [0; 0; 0];   % torque contribution (added after the strip loop)
F_span_full     = [0; 0; 0];   % full span force, for intermediates (zero when disabled)
spanTorqueAtten = 1;           % span-torque reduced-frequency attenuation factor
r_spanCoP_body  = [0; 0; 0];   % span-CoP application point relative to CoM, body frame

if enableSpanForce
    % Whole-seed geometric centre in body coords (area-weighted; body y = 0).
    seedGeoCenter = [xGeoCenterChord; 0; zGeoCenterSpan];

    % Bulk velocity at the geometric centre. With the transport term (default):
    %   v_gc^body = R' * v_com^inertial + omega^body x (geoCentre - CoM)
    % Without it, only the CoM translational velocity drives the span force, so
    % the seed's own rotational sweep cannot feed back into this translational
    % drag (see the toggle notes above).
    v_com_body = R.' * v;
    if useSpanGeomVelocity
        v_gc_body = v_com_body + cross(omega, seedGeoCenter - mp.c);
    else
        v_gc_body = v_com_body;      % translation only -- omit the rotational sweep
    end
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
    %
    % Application point is switchable: the migrating span CoP (default) or the
    % fixed geometric centre. See the enableSpanCOPMigration notes above.
    if useSpanCOPMigration
        zSpanCoP  = computeStripCoP(l_cp_frac_span, totalSpan, zGeoCenterSpan);  % migrating span CoP (body z)
        r_spanArm = [xGeoCenterChord; 0; zSpanCoP] - mp.c;   % span-CoP arm, relative to CoM
    else
        r_spanArm = seedGeoCenter - mp.c;   % fixed geometric-centre arm (no CoP migration)
    end
    r_spanCoP_body = r_spanArm;   % expose the span-CoP application point (rel. to CoM)

    % Reduced-frequency attenuation (quasi-steady validity). The span torque is a
    % lumped quasi-steady CoP moment, valid only when the span-plane flow direction
    % changes slowly relative to the convective time S/v_ip. During autorotation the
    % descent velocity resolved into the spinning body frame makes beta -- and this
    % torque -- oscillate at the NORMAL-axis spin omega_y, parametrically
    % destabilising roll; the reduced frequency k = |omega_y|*S/(2*v_ip) captures
    % that. Keying on omega_y (not total |omega|) leaves the tumbling mode (an
    % omega_z mode, small omega_y) unattenuated.
    if useSpanTorqueAttenuation
        v_ip_span       = hypot(vSpan_gc, vNormal_gc);
        k_reduced       = abs(omega_y) * totalSpan / (2 * max(v_ip_span, eps));
        spanTorqueAtten = 1 / (1 + (k_reduced / spanCoeffs.k0_spanTorque)^2);
    else
        spanTorqueAtten = 1;   % raw, unattenuated span torque
    end

    % C_span_torque scales the torque only, leaving F_span_apply untouched, so
    % the force and torque strengths are independently tunable. enableSpanTorque
    % gates the moment entirely (force still applied) for on/off debugging.
    if useSpanTorque
        tau_span = spanTorqueAtten * spanCoeffs.C_span_torque * cross(r_spanArm, F_span_full);
    else
        tau_span = [0; 0; 0];   % span FORCE kept, span TORQUE dropped
    end
    %tau_span = spanTorqueAtten * spanCoeffs.C_span_torque * cross(r_spanArm, F_span_apply);
end

for i = 1 : numStrips

    % --- strip local frame (columns of strip->body rotation) and 3D centre
    cDir = chordDir(:, i);   nDir = normalDir(:, i);   sDir = spanDir(:, i);
    pGeo = [xGeo(i); yGeo(i); zGeo(i)];          % strip geometric centre, body coords

    % --- strip velocity, projected into the strip's LOCAL (chord, normal) plane.
    % vChord/vNormal are the components the 2D section sees; the local-span
    % component is discarded (as the planar model discards body-z). For a flat
    % plate cDir=[1;0;0], nDir=[0;1;0] -> vChord=vStrip(1), vNormal=vStrip(2).
    vStrip  = stripVel(:, i);
    vChord  = cDir.' * vStrip;
    vNormal = nDir.' * vStrip;

    % --- angle of attack in the local plane (same law, projected components)
    alpha = computeAngleOfAttack([vChord; vNormal; 0]);

    % --- aerodynamic coefficients at this angle of attack
    coeffs = computeAeroCoeffs(alpha, aeroParams);

    % --- spin about the strip's LOCAL span axis drives the 2D rotational lift/damp
    omega_localSpan = sDir.' * omega;

    % --- forces in the strip's LOCAL frame ([chord; normal; 0]), then rotate to
    % body via the strip frame R_strip = [cDir nDir sDir]. Per-strip liftMult/
    % dragMult scale translational lift and drag only (not the rotational lift).
    [F_transl_L, F_rotLift_L] = computeStripForces(vChord, vNormal, chord(i), dz(i), ...
                                               omega_localSpan, coeffs, rhoFluid, ...
                                               liftMult(i), dragMult(i));
    F_transl  = cDir*F_transl_L(1)  + nDir*F_transl_L(2)  + sDir*F_transl_L(3);
    F_rotLift = cDir*F_rotLift_L(1) + nDir*F_rotLift_L(2) + sDir*F_rotLift_L(3);
    dF = F_transl + F_rotLift;

    % --- centre of pressure: offset l_cp_frac*chord along the LOCAL chord axis
    % from the strip centre. (xcp is the equivalent chordwise coordinate, kept
    % for the intermediates.)
    cpOffset    = coeffs.l_cp_frac * chord(i);
    xcp         = xGeo(i) + cpOffset;                 % = computeStripCoP(...)
    r_cp        = (pGeo + cpOffset*cDir) - mp.c;       % CoP relative to CoM
    r_geoCenter = pGeo - mp.c;                         % geometric centre rel. to CoM

    % --- torques about the CoM: each force at its own point
    tau_transl  = cross(r_cp,        F_transl);   % 2D Tt   (at CoP)
    tau_rotLift = cross(r_geoCenter, F_rotLift);  % 2D Tcr  (at geometric centre / mid-chord)

    % --- extra torque: nonlinear spin damping about the LOCAL span axis (2D Tr).
    % CoM offset measured along the local chord from the strip centre.
    comOffsetFromGeoCenter = cDir.' * (mp.c - pGeo);
    Tr        = stripSpinDamping(chord(i), dz(i), omega_localSpan, comOffsetFromGeoCenter, coeffs.CD_rot, rhoFluid);
    tau_spin  = Tr * sDir;                        % damping about the local span axis

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
%   - spin-damping torques Tx (chordwise/roll, gated + scaled) and Ty (normal)
%   - the optional spanwise-flow force (body-z only) and its torque
% Tx is off by default (see enableTxDamping); when on it is scaled by C_Tx.
if useTxDamping
    Tx_applied = C_Tx_const * Tx;
else
    Tx_applied = 0;
end
tau_body    = tau_body    + [Tx_applied; Ty; 0] + tau_span;
F_aero_body = F_aero_body + F_span_apply;

% =========================================================================
% 4. SHARED RIGID-BODY CORE  ->  state derivative
% =========================================================================
% The shape-specific aerodynamics are done above (F_aero_body, tau_body). Hand
% the net body-frame force/torque + mass properties to the shape-agnostic 6-DOF
% core, which rotates the force to inertial, adds gravity, solves the
% translational (with added mass) and modified-Euler rotational dynamics, and
% integrates the quaternion. The SAME core is used by the 3D-shape force model.
[xdot, core] = rigidBody6DOF(x, mp, F_aero_body, tau_body, gMag);

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
    intermediates.F_aero_inertial = core.F_aero_inertial;  % 3x1 total aero force, inertial
    intermediates.F_total_inertial= core.F_total_inertial; % 3x1 aero + gravity, inertial
    intermediates.Tx_spanSpin     = Tx;               % computed chordwise (roll) spin-damping torque (N*m)
    intermediates.Tx_applied      = Tx_applied;       % Tx actually added (0 if enableTxDamping false)
    intermediates.Ty_normalSpin   = Ty;               % whole-seed normal-axis spin-damping torque (N*m)
    intermediates.spanTorqueAtten = spanTorqueAtten;  % span-torque reduced-frequency attenuation factor
    intermediates.F_span_full     = F_span_full;      % 3x1 full span-flow force, body (y-component discarded from sum)
    intermediates.F_span_apply    = F_span_apply;     % 3x1 span-flow force actually added (body-z only)
    intermediates.r_spanCoP_body  = r_spanCoP_body;   % 3x1 span-CoP application point relative to CoM, body frame
    intermediates.tau_span        = tau_span;         % 3x1 span-flow torque, body (uses full force + migrating span CoP)
    intermediates.tau_body        = tau_body;         % 3x1 total torque, body
    intermediates.a_inertial      = core.a_inertial;  % 3x1 linear acceleration
    intermediates.alpha_body      = core.alpha_body;  % 3x1 angular acceleration
end

end
