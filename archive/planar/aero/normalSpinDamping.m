function Ty = normalSpinDamping(R, omega_y, CD0, C_fy, rhoFluid)
% NORMALSPINDAMPING  Spin-damping torque about the PLATE-NORMAL (body y) axis,
% for the whole seed.
%
% Placeholder physics: in-plane (normal-axis) spin sweeps the plate edge-on to
% the flow along both the chord and span -- the same motion getAddedMass.m
% treats as ~zero added mass (m1 = 0) -- so this is NOT a rederivation of the
% stripSpinDamping/spanSpinDamping chord/span-end quartic integral; there is
% no analogous clean closed form for this axis yet. Instead it borrows the
% SAME -1/128*rho*(length)^4 scaling (for magnitude consistency with the
% other two damping torques) and the existing zero-incidence drag coefficient
% CD0 (a real, measured quantity from computeAeroCoeffs, rather than a
% fabricated constant), multiplied by a free tuning factor C_fy so normal-
% axis spin can be made much easier to excite than spin about x or z, as
% expected for a near-flat plate.
%
% INPUTS
%   R        : characteristic radius (m). Currently: half the seed's total span.
%   omega_y  : plate-normal (body y) angular velocity (rad/s).
%   CD0      : zero-incidence/parasitic drag coefficient (from computeAeroCoeffs).
%   C_fy     : tuning factor for this damping mode (dimensionless; aero.C_fy,
%              default 1.0 -- decrease to make normal-axis spin freer).
%   rhoFluid : fluid density (kg/m^3).
%
% OUTPUT
%   Ty : spin-damping torque about the plate-normal (body y) axis (N*m). Opposes spin.

    Ty = -1/128 * rhoFluid * R^4 * CD0 * C_fy * omega_y * abs(omega_y);
end
