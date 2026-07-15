function [F_span_full, l_cp_frac_span] = computeSpanForce(vSpan, vNormal, totalSpan, chordMean, C_span, coeffs, rhoFluid)
% COMPUTESPANFORCE  Whole-seed aerodynamic force from SPANWISE flow, body frame.
%
% The chordwise strips (computeStripForces) model flow in the chord-normal
% (body x-y) plane and produce ZERO spanwise (body z) force by construction, so
% a seed sliding sideways feels no aerodynamic resistance. This function fills
% that gap with a single whole-seed "span strip": the exact same quasi-steady
% lift/drag model, but applied in the span-normal (body z-y) plane, with the
% SPAN playing the "chord" role and the mean chord playing the strip-width role
% (so span * chordMean = total wing area, the area presented to spanwise flow).
%
% It returns the FULL force (including a normal / body-y component). The caller
% is expected to DISCARD the body-y component from the net force sum -- the
% strips already model the normal-direction force, and adding it again here
% would double-count it. See the caller (seed6DOFODE) and the notes below for
% exactly what is discarded and why the full force is nonetheless returned.
%
% Direction convention (mirrors computeStripForces, with the chord axis x
% replaced by the span axis z): in-plane velocity (vSpan, vNormal) gives
% lift ~ (vNormal, -vSpan) and drag ~ (vSpan, vNormal), placed into the body
% (z, y) components; the body-x (chord) component is identically zero.
%
% INPUTS
%   vSpan     : spanwise (body z) velocity of the seed geometric centre (m/s).
%   vNormal   : plate-normal (body y) velocity of the seed geometric centre (m/s).
%   totalSpan : total spanwise extent of the seed (m). Plays the "chord" role.
%   chordMean : mean aerodynamic chord (m). Plays the strip-width ("dz") role.
%   C_span    : spanwise-flow force tuning factor (dimensionless; aero.C_span,
%               default 1.0 -- scale the whole span force up/down while tuning).
%   coeffs    : struct from computeAeroCoeffs evaluated at the SPAN-plane angle
%               of attack beta = atan2(vNormal, vSpan). Uses .CT, .CD, .l_cp_frac.
%   rhoFluid  : fluid density (kg/m^3).
%
% OUTPUTS
%   F_span_full    : 3x1 span-flow force, body frame (N), = [0; Fy; Fz].
%                    The body-x (chord) component is 0. The body-y (normal)
%                    component IS populated here but is meant to be discarded
%                    from the net force by the caller (see above).
%   l_cp_frac_span : centre-of-pressure offset as a FRACTION of the SPAN, from
%                    the span geometric centre (the span-plane analogue of the
%                    strip l_cp_frac). The caller turns this into a migrating
%                    span-CoP position for the torque arm.

    % In-plane speed in the span-normal plane (spanwise + normal only).
    v_ip = sqrt(vSpan^2 + vNormal^2);

    % Magnitude factors (each carries one power of speed). Span plays the
    % "chord" role and chordMean the "dz" role, so span*chordMean = wing area.
    % C_span scales the whole thing for tuning.
    Lt_0 = 0.5 * rhoFluid * totalSpan * chordMean * coeffs.CT * v_ip * C_span;   % lift  factor
    D_0  = -0.5 * rhoFluid * totalSpan * chordMean * coeffs.CD * v_ip * C_span;  % drag  factor (opposes v)

    % Assemble in the body frame. This is computeStripForces' F_transl with the
    % chord axis (x) swapped for the span axis (z):
    %   lift ~ (vNormal, -vSpan) ,  drag ~ (vSpan, vNormal)  in the (z, y) plane.
    F_span_full = [ 0 ;                              % body x (chord): span flow exerts no chordwise force
                   -Lt_0*vSpan  + D_0*vNormal ;      % body y (normal): populated, but CALLER DISCARDS from force sum
                    Lt_0*vNormal + D_0*vSpan ];       % body z (span):  the previously-missing spanwise force

    % CoP migrates along the SPAN (body z) with beta, exactly as the strip CoP
    % migrates along the chord with alpha. Returned as a fraction of span.
    l_cp_frac_span = coeffs.l_cp_frac;
end
