function xcp = computeStripCoP(l_cp_frac, chord, xMid)
% COMPUTESTRIPCOP  Chordwise center-of-pressure position for one strip (body x).
%
% The center of pressure migrates along the chord with angle of attack. Its
% offset from the strip mid-chord is l_cp_frac * chord (l_cp_frac from
% computeAeroCoeffs). Passing l_cp_frac = 0 recovers the simplest "center of
% section" model (CoP at mid-chord).
%
% INPUTS
%   l_cp_frac : center-of-pressure offset as a fraction of chord, from mid-chord.
%   chord     : local strip chord length (m).
%   xMid      : chordwise position of the strip mid-chord in the body frame (m)
%               (0 for the symmetric seed).
%
% OUTPUT
%   xcp : chordwise CoP position in the body frame (m).

    xcp = xMid + l_cp_frac * chord;
end
