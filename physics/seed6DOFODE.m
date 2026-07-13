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
    [F_transl, F_rotLift] = computeStripForces(vChord, vNormal, chord(i), dz(i), ...
                                               omega_z, coeffs, rhoFluid);
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
    intermediates.tau_body        = tau_body;         % 3x1 total torque, body
    intermediates.a_inertial      = a_inertial;       % 3x1 linear acceleration
    intermediates.alpha_body      = alpha_body;       % 3x1 angular acceleration
end

end
