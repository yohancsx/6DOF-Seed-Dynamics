function [F_transl, F_rotLift] = computeStripForces(vChord, vNormal, chord, dz, omega_z, coeffs, rhoFluid, liftMult, dragMult)
% COMPUTESTRIPFORCES  Aerodynamic forces on one strip, body frame.
%
% Reproduces the 2D quasi-steady force model (minimal_imp lines 174-185) per
% strip, splitting the force into:
%   F_transl  : translational lift + drag (from CT, CD)
%   F_rotLift : rotational (added-circulation) lift (from CR and spanwise spin)
% They are returned separately because their TORQUES act at different chordwise
% points: translational at the CoP, rotational lift at mid-chord. Both lifts
% point perpendicular to the in-plane velocity, in the same sense; drag opposes
% the velocity.
%
% Direction convention (from the 2D code): with in-plane velocity (vx, vy),
% lift ~ (vy, -vx) and drag ~ (vx, vy). Multiplying a "magnitude factor" (which
% already carries one power of speed) by the velocity supplies the second speed
% power and the direction.
%
% INPUTS
%   vChord   : chordwise (body x) velocity at the strip (m/s).
%   vNormal  : plate-normal (body y) velocity at the strip (m/s).
%   chord    : strip chord length (m).
%   dz       : strip spanwise width (m).
%   omega_z  : spanwise (body z) angular velocity component (rad/s) -- the only
%              spin that drives the 2D rotational lift.
%   coeffs   : struct from computeAeroCoeffs (uses .CT, .CD, .CR).
%   rhoFluid : fluid density (kg/m^3).
%   liftMult : (optional) static multiplier on this strip's TRANSLATIONAL lift,
%              a per-strip tuning/morphology knob. Default 1 (no change). Does
%              NOT scale the rotational lift.
%   dragMult : (optional) static multiplier on this strip's drag. Default 1.
%
% OUTPUTS
%   F_transl  : 3x1 translational lift + drag, body frame (N). z-component 0.
%   F_rotLift : 3x1 rotational lift,           body frame (N). z-component 0.

    % Per-strip lift/drag multipliers default to 1 (no scaling) when omitted.
    if nargin < 8 || isempty(liftMult); liftMult = 1; end
    if nargin < 9 || isempty(dragMult); dragMult = 1; end

    % In-plane speed (spanwise component already discarded upstream).
    v_ip = sqrt(vChord^2 + vNormal^2);

    % Magnitude factors (each carries one power of speed; 0.5*rho*area = 0.5*rho*chord*dz).
    % liftMult/dragMult scale the translational lift and drag only; the rotational
    % lift factor Lr_0 is deliberately left unscaled.
    Lt_0 = liftMult * 0.5 * rhoFluid * chord * dz * coeffs.CT * v_ip;   % translational lift [2D L_0, CT part]
    D_0  = -dragMult * 0.5 * rhoFluid * chord * dz * coeffs.CD * v_ip;  % drag (opposes v)   [2D D_0]
    Lr_0 = -0.5 * rhoFluid * chord^2 * dz * coeffs.CR * omega_z;        % rotational lift    [2D L_0, CR part]

    % Assemble in the body frame. Lift ~ (vNormal, -vChord); drag ~ (vChord, vNormal).
    F_transl  = [ Lt_0*vNormal + D_0*vChord ;
                 -Lt_0*vChord  + D_0*vNormal ;
                  0 ];

    F_rotLift = [ Lr_0*vNormal ;
                 -Lr_0*vChord  ;
                  0 ];
end
