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
% PHASE 4 -- the span force went too (measured inert, wrong shape), and two
% physically-grounded terms came in, each behind a switch that setupSeedShape3D
% turns ON by default:
%   * enableEdgeDrag      -- pure drag on a seed sliding along its span, which
%     the strips cannot see at all (computeEdgeDrag). Section 3b.
%   * enableAddedMassRate -- the Adot*v term of the translational EOM: spinning
%     changes the entrained fluid's momentum even at constant speed. Lives in the
%     shared core, rigidBody6DOF; off there by default so planar stays frozen.
%
% FLAT EQUIVALENCE: with the planar extras disabled (enableSpanForce,
% enableTxDamping, enableNormalSpinDamping false on the planar side) and the two
% phase-4 terms disabled here (enableEdgeDrag, enableAddedMassRate false), this
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
% 3a. WHOLE-SEED GEOMETRY AGGREGATES  (used by the edge drag in 3b)
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
yGeoCenterNorm  = sum(chord .* dz .* yGeo) / wingArea; % area-weighted out-of-plane
                                                        % centre (m); 0 unless curved

% =========================================================================
% 3b. EDGE CROSSFLOW DRAG  (enableEdgeDrag, DEFAULT TRUE for shape3d)
% =========================================================================
% The chordwise strips produce exactly zero force along the span, so a seed
% sliding sideways would feel no aerodynamic resistance at all. Physically it
% presents its wingtip -- a small bluff section, thickness by chord. See
% computeEdgeDrag for the law and why the area centroid is the exact application
% point. Pure drag opposing the spanwise velocity; no lift, no migrating CoP.
%
% The velocity is sampled at the application point by the transport theorem,
% v = R'*v_CoM + omega x (centroid - CoM). That is simply the correct velocity at
% that point: the retired span force made this optional (enableSpanGeomVelocity)
% and defaulted it OFF, which was the unphysical choice.
%
% RETIRED HERE (phase 4): computeSpanForce and its two switches, enableSpanForce
% and enableSpanGeomVelocity. Phase 3 measured the span force inert (<= 7e-4
% relative on every benchmark metric), and its law had the wrong shape -- smallest
% for a pure spanwise slide and largest where the strips already supply the force.
% Its torque apparatus had already gone in the Sep-2026 audit. The planar model
% keeps it unchanged. This term is expected to be inert too (~4x smaller again);
% it exists to remove a known-unphysical zero, not to move results.
useEdgeDrag = isfield(seedParams, 'enableEdgeDrag') && seedParams.enableEdgeDrag;

% Contributions default to zero so the post-loop accumulation is unconditional.
F_edge         = [0; 0; 0];    % edge-drag force (added after the strip loop), body
tau_edge       = [0; 0; 0];    % its moment about the CoM, body
r_edgeArm_body = [0; 0; 0];    % application point relative to the CoM, body
vSpan_edge     = 0;            % spanwise velocity the drag responded to (m/s)

if useEdgeDrag
    areaCentroid   = [xGeoCenterChord; yGeoCenterNorm; zGeoCenterSpan];
    r_edgeArm_body = areaCentroid - mp.c;
    v_centroid     = R.' * v + cross(omega, r_edgeArm_body);   % transport theorem
    A_edge         = seedParams.seedThickness * chordMean;     % frontal area (m^2)
    edgeCoeffs     = computeAeroCoeffs(0, aeroParams);         % C_d_edge is alpha-free
    [F_edge, tau_edge, vSpan_edge] = computeEdgeDrag(v_centroid, r_edgeArm_body, ...
        [0; 0; 1], A_edge, edgeCoeffs.C_d_edge, rhoFluid);
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
% strip-accumulation loop above: the edge crossflow drag and its moment.
%
% The whole-seed spin-damping torques Tx and Ty that used to be added here are
% gone (see the header): Tx duplicated the strip loop's own roll moment, and Ty
% was dimensionally a force. The span force is gone too, replaced by edge drag.
tau_body    = tau_body    + tau_edge;
F_aero_body = F_aero_body + F_edge;

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
    % r_spanCoP_body, and the span-force fields F_span_full / F_span_apply /
    % tau_span, are GONE -- the terms they reported were removed (see the header).
    % The planar seed6DOFODE still exposes them, and their only consumer,
    % testing/planar/runSpanForceComparison.m, is a planar-only diagnostic.
    intermediates.F_edge          = F_edge;           % 3x1 edge crossflow drag, body (0 if off)
    intermediates.tau_edge        = tau_edge;         % 3x1 its moment about the CoM, body
    intermediates.r_edgeArm_body  = r_edgeArm_body;   % 3x1 area centroid relative to CoM, body
    intermediates.vSpan_edge      = vSpan_edge;       % spanwise velocity at the centroid (m/s)
    intermediates.F_addedMassRate = core.F_addedMassRate; % 3x1 Adot*v, inertial (0 if off)
    intermediates.tau_body        = tau_body;         % 3x1 total torque, body
    intermediates.a_inertial      = core.a_inertial;  % 3x1 linear acceleration
    intermediates.alpha_body      = core.alpha_body;  % 3x1 angular acceleration
end

end
