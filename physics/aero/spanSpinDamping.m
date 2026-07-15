function Tx = spanSpinDamping(totalSpan, chord, omega_x, comSpanOffset, CD_rot, rhoFluid)
% SPANSPINDAMPING  Spin-damping torque about the CHORDWISE (body x) axis, for
% the WHOLE seed.
%
% Direct analogue of stripSpinDamping (2D Tr, minimal_imp line 195), with the
% integrated axis and the multiplicative-width axis swapped. Rotation about
% body x induces a velocity ~ omega_x*(z-zcm) that varies with SPANWISE
% position, not chordwise position (the opposite of the omega_z case), so the
% same quartic chord-end integral is redone here over the SPAN, using the
% seed's total span in the "l" role and the chord as the simple multiplicative
% width (the role dz played in stripSpinDamping).
%
% IMPORTANT -- evaluate ONCE per timestep, not per strip:
% stripSpinDamping is called once per strip because its chordwise integral is
% identical at every spanwise position (for a non-tapered wing), so summing
% each strip's dz-weighted contribution exactly reconstructs the full
% chord-integral x total-span product. Here it is the opposite: span IS the
% axis the strip loop already discretises, so this function performs that
% span integral itself in closed form. Calling it once per strip with the
% same totalSpan would replay the whole-span integral redundantly and wildly
% overcount the torque -- call it exactly once, outside the strip loop.
%
% INPUTS
%   totalSpan     : total spanwise extent of the seed (m). The "l" of the
%                   quartic chord-end integral, now taken over the span.
%   chord         : representative chord (m) -- the simple multiplicative
%                   width (dz's role in stripSpinDamping). Exact for a
%                   non-tapered wing (constant chord); for a tapered wing use
%                   the mean aerodynamic chord, wingArea / totalSpan.
%   omega_x       : chordwise (body x) angular velocity (rad/s).
%   comSpanOffset : spanwise offset of the CoM from the seed's OVERALL
%                   spanwise geometric centre (m) -- the "lcm" analogue, now
%                   measured along the span instead of the chord.
%   CD_rot        : rotational damping coefficient (= CD2 from computeAeroCoeffs).
%   rhoFluid      : fluid density (kg/m^3).
%
% OUTPUT
%   Tx : spin-damping torque about the chordwise (body x) axis (N*m). Opposes spin.

    ratio     = 2 * comSpanOffset / totalSpan;         % 2*lcm/l analogue, span version
    endFactor = (ratio + 1)^4 + (ratio - 1)^4;         % span-end integration limits

    Tx = -1/128 * rhoFluid * totalSpan^4 * chord * CD_rot ...
         * omega_x * abs(omega_x) * endFactor;
end
