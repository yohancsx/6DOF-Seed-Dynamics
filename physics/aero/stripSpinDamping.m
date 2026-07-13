function Tr = stripSpinDamping(chord, dz, omega_z, xcm_local, CD_rot, rhoFluid)
% STRIPSPINDAMPING  Spin-damping torque about the spanwise axis for one strip.
%
% Generalises the 2D Tr term (minimal_imp line 195):
%   Tr = -1/128 * rho * l^4 * CD2 * omega|omega| * ((2lcm/l+1)^4 + (2lcm/l-1)^4)
% As the plate spins, each chordwise element feels drag ~ (omega*r)^2; integrating
% r x drag along the chord about the offset CoM gives this nonlinear (omega|omega|)
% damping. The two (2*lcm/l +/- 1)^4 terms are the chord-end integration limits.
% Per strip, l -> chord, add the strip width dz, and use the CoM's chordwise
% offset from THIS strip's mid-chord as the lever offset.
%
% This is separate physics from the strip lift/drag, which are evaluated at
% mid-chord where the spanwise spin contributes no velocity.
%
% INPUTS
%   chord     : strip chord length (m).
%   dz        : strip spanwise width (m).
%   omega_z   : spanwise (body z) angular velocity (rad/s).
%   xcm_local : chordwise offset of the CoM from the strip mid-chord (m)
%               (the per-strip analogue of the 2D lcm; 0 for the symmetric seed).
%   CD_rot    : rotational damping coefficient (= CD2 from computeAeroCoeffs).
%   rhoFluid  : fluid density (kg/m^3).
%
% OUTPUT
%   Tr : spin-damping torque about the spanwise (body z) axis (N·m). Opposes spin.

    ratio     = 2 * xcm_local / chord;                 % 2*lcm/l analogue
    endFactor = (ratio + 1)^4 + (ratio - 1)^4;         % chord-end integration limits

    Tr = -1/128 * rhoFluid * chord^4 * dz * CD_rot ...
         * omega_z * abs(omega_z) * endFactor;
end
