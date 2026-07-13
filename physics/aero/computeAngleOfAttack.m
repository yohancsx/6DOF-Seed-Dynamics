function alpha = computeAngleOfAttack(vStripBody)
% COMPUTEANGLEOFATTACK  Local angle of attack for one strip.
%
% Uses the strip's in-plane velocity components (chordwise x, normal y) in the
% body frame. Matches the 2D model's definition alpha = atan2(vyp, vxp) and the
% convention expected by computeAeroCoeffs (velocity-based, NOT wind-based). The
% spanwise (z) component is ignored -- strip theory treats each strip as 2D.
%
% INPUT
%   vStripBody : 3x1 body-frame velocity at the strip [chordwise; normal; spanwise].
%
% OUTPUT
%   alpha : angle of attack (rad), measured from the chord (body x) to the
%           in-plane velocity.

    vChord  = vStripBody(1);   % body x
    vNormal = vStripBody(2);   % body y
    alpha   = atan2(vNormal, vChord);
end
