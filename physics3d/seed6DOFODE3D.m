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
% SCOPE -- WHAT THIS MODEL DELIBERATELY CONTAINS:
%   Andersen-Pesavento-Wang sectional coefficients  +  strip theory  +  rigid-body
%   first principles  +  added mass.  Nothing else.
%
% A Sep-2026 audit removed the terms that had been invented to make a FLAT strip
% decomposition fake out-of-plane behaviour. They are gone from this RHS because a
% seed with twist and curvature produces that behaviour from real geometry:
%   * Tx (spanSpinDamping) -- provably a DOUBLE COUNT. A pure-roll state isolates
%     the strips' own roll moment and it converges on Tx (ratio 1.0001 at 192
%     strips); both reduce to -rho*c*CD2*w|w|*S^4/64. Span IS the axis the strip
%     loop discretises, so the closed-form span integral re-adds what the strips
%     already computed. Roll damping was 2x reality.
%   * Ty (normalSpinDamping) -- dimensionally a FORCE (rho*R^4*w^2 = N, not N*m;
%     Tr and Tx each carry a width that Ty lacked) and 0.1% of the torque budget.
%   * The span TORQUE, its migrating span-CoP, and the reduced-frequency
%     attenuation -- 25% of the torque budget generated from F_span_full(2), a
%     force this code explicitly EXCLUDES from the net force as a double count.
%     A moment with no corresponding force matches no pressure distribution.
%   * The tuning constants those terms carried: C_span_torque, k0_spanTorque,
%     C_Tx, C_fy. They still exist in the shared computeAeroCoeffs (the frozen
%     planar model reads them); this RHS simply never asks for them.
% The planar seed6DOFODE keeps all of the above, unchanged, as the frozen
% reference model. See the README "Model cleanup" roadmap section.
%
% The span FORCE survives (enableSpanForce), but now contributes its moment the
% honest way -- cross(arm, force_actually_applied) at the seed geometric centre.
% setupSeedShape3D switches it off for a curved seed, whose tilted strips supply
% that force from geometry.
%
% FLAT EQUIVALENCE: with the planar extras disabled (enableSpanForce,
% enableTxDamping, enableNormalSpinDamping all false on the planar side) this
% still reduces EXACTLY to seed6DOFODE on a flat seed -- see
% testing/shape3d/testShape3DFlatEquivalence.m.
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
% 3a. WHOLE-SEED GEOMETRY AGGREGATES  (used by the span force in 3b)
% =========================================================================
% NOTE: the whole-seed spin-damping torques Tx (spanSpinDamping) and Ty
% (normalSpinDamping) that used to live here have been REMOVED -- see the header.
% Tx duplicated the roll moment the strip loop already produces, and Ty was
% dimensionally a force. Spin damping about the strip span axis is still present,
% per strip, via stripSpinDamping (2D Tr) inside the loop below: that one is real
% sub-strip physics, because the strip loop does not resolve the CHORD.
totalSpan = sum(dz);               % total spanwise extent of the seed (m)
wingArea  = sum(chord .* dz);      % total wing area (m^2)
chordMean = wingArea / totalSpan;  % mean aerodynamic chord (m); exact for a
                                    % non-tapered wing, representative otherwise
zGeoCenterSpan = sum(chord .* dz .* zGeo) / wingArea;  % area-weighted spanwise
                                                        % geometric centre (m)
xGeoCenterChord = sum(chord .* dz .* xGeo) / wingArea; % area-weighted chordwise
                                                        % geometric centre (m)

% =========================================================================
% 3b. WHOLE-SEED SPANWISE-FLOW FORCE  (optional; fills the unmodeled body-z force)
% =========================================================================
% The chordwise strips produce zero body-z (spanwise) force, so a sideways-
% sliding seed feels no aerodynamic resistance. computeSpanForce supplies a
% single whole-seed force from flow in the span-normal (z-y) plane.
%
% Two switches (both behave as in the planar model when the field is absent):
%   .enableSpanForce        - apply the span force at all (DEFAULT TRUE).
%                             setupSeedShape3D sets it FALSE for a CURVED seed,
%                             whose tilted strips make this force from geometry.
%   .enableSpanGeomVelocity - include the rotational omega x r transport when
%                             sampling the velocity at the geometric centre
%                             (DEFAULT FALSE, matching the planar model).
%     The omega x r term feeds a SPANWISE velocity only when the
%     geometric-centre-to-CoM arm has a CHORDWISE (or out-of-plane) component,
%     since (omega x r)_z = omega_x*r_y - omega_y*r_x. For a purely spanwise CoM
%     offset it contributes nothing, so toggling this isolates how much the
%     rotational sweep (vs. the descent velocity resolved in the spinning body
%     frame) drives the span force.
%
% REMOVED (Sep-2026 audit) -- the span force's whole torque apparatus:
%   .enableSpanCOPMigration, .enableSpanTorque, .enableSpanTorqueAttenuation and
%   the constants C_span_torque / k0_spanTorque. The old span torque crossed a
%   migrating span-CoP arm with the FULL span force, whose body-y component is
%   excluded from the net force as a double count -- measured, 100% of that
%   torque (25% of the seed's whole moment budget) came from a force never
%   applied. The moment below is now just cross(arm, force actually applied).
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

% Contributions default to zero so the post-loop accumulation is unconditional.
F_span_apply    = [0; 0; 0];   % force  contribution (added after the strip loop)
tau_span        = [0; 0; 0];   % torque contribution (added after the strip loop)
F_span_full     = [0; 0; 0];   % full span force, for intermediates (zero when disabled)
r_spanApply_body = [0; 0; 0];  % span-force application point relative to CoM, body frame

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

    % Full span force [0; Fy; Fz]. (computeSpanForce also returns a span-CoP
    % fraction; it is no longer used -- see the CoP-migration note in the header.)
    F_span_full = computeSpanForce(vSpan_gc, vNormal_gc, ...
        totalSpan, chordMean, spanCoeffs.C_span, spanCoeffs, rhoFluid);

    % --- FORCE: keep ONLY the body-z (spanwise) component ------------------
    % Explicitly DISCARD body-x (already 0) and body-y (normal). The strips
    % already model the normal-direction force; re-adding F_span_full(2) would
    % double-count it. Only the body-z direction was previously unmodeled.
    F_span_apply = [0; 0; F_span_full(3)];

    % --- TORQUE: the moment of the force that is ACTUALLY APPLIED ----------
    % A force acting at a point offset from the CoM produces a moment; that is the
    % whole of it. Applied at the seed's area-weighted geometric centre, which is
    % where a uniformly-distributed spanwise load acts.
    %
    % This replaces the old "Option B" torque, which crossed a migrating span-CoP
    % arm with the FULL span force -- including the body-y component deliberately
    % excluded from the net force two lines above. Measured over a settled
    % autorotation, 100% of that torque came from the excluded component (by a
    % factor of 640), and it was 25% of the seed's entire moment budget. A moment
    % with no corresponding force is consistent with no pressure distribution, so
    % linear and angular momentum were being fed different loads. The stated
    % justification -- that strip theory cannot produce a roll moment from
    % span-offset loading -- is also false here: each strip's normal force acts at
    % its own zGeo, so cross(r_cp, F) already carries an x-component whenever
    % zGeo ~= c_z.
    r_spanApply_body = seedGeoCenter - mp.c;
    tau_span         = cross(r_spanApply_body, F_span_apply);
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

% Add the whole-seed contribution (section 3b) -- a single contribution for the
% whole seed, not a per-strip one, so it is added here rather than inside the
% strip-accumulation loop above: the optional spanwise-flow force (body-z only)
% and the moment it makes at its application point.
%
% The whole-seed spin-damping torques Tx and Ty that used to be added here are
% gone (see the header): Tx duplicated the strip loop's own roll moment, and Ty
% was dimensionally a force.
tau_body    = tau_body    + tau_span;
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
    % NOTE: Tx_spanSpin / Tx_applied / Ty_normalSpin / spanTorqueAtten /
    % r_spanCoP_body are GONE -- the terms they reported were removed (see the
    % header). The planar seed6DOFODE still exposes them, and the only consumer,
    % testing/planar/runSpanForceComparison.m, is a planar-only diagnostic.
    intermediates.F_span_full     = F_span_full;      % 3x1 full span-flow force, body (y-component discarded from sum)
    intermediates.F_span_apply    = F_span_apply;     % 3x1 span-flow force actually added (body-z only)
    intermediates.r_spanApply_body= r_spanApply_body; % 3x1 span-force application point relative to CoM, body frame
    intermediates.tau_span        = tau_span;         % 3x1 span-flow torque, body = cross(arm, F_span_apply)
    intermediates.tau_body        = tau_body;         % 3x1 total torque, body
    intermediates.a_inertial      = core.a_inertial;  % 3x1 linear acceleration
    intermediates.alpha_body      = core.alpha_body;  % 3x1 angular acceleration
end

end
