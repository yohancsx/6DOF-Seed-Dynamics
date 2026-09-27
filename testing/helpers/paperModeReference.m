function mode = paperModeReference(xa, yb)
% PAPERMODEREFERENCE  The published mode at a point of Hou et al.'s Fig. 2a.
%
% A HAND-DIGITISED reading of their phase map, so a simulated grid can be scored
% against it cell by cell. Returns 'ST' | 'CH' | 'AR' | 'FA' for a point
% (x_c/a, y_c/b).
%
% HONESTY ABOUT WHAT THIS IS. The boundaries below were read off the published
% figure by eye, at figure resolution, and the regions there are drawn as smooth
% shaded bands rather than sharp lines. Treat an agreement score from this as
% meaningful to a few per cent at best, and treat any cell within ~0.02 of a
% boundary as unresolved. The CONTINUOUS vs SEGMENTED split inside ST (their
% dashed line) is NOT encoded: at this resolution I could not place it reliably
% enough to score against, so ST is returned as one region.
%
% The boundaries, in the order they are tested:
%   FA   x >= 0.28  and  y <= 0.115          (their purple block, lower right)
%   AR   y >= max(0.115, 0.205 - 0.40*x)     (yellow; the lower edge rises to the
%                                             left, reaching ~0.19 at x ~ 0.05)
%   ST   y <= 0.145 - 0.50*x                 (red; lower-left triangle, reaching
%                                             y ~ 0.145 at x = 0 and y ~ 0 at x ~ 0.29)
%   CH   everything else                     (blue; the diagonal band between)
%
% INPUTS   xa = x_c/a (spanwise), yb = y_c/b (chordwise). Scalars or arrays.
% OUTPUT   mode : char for scalar input, cellstr for array input.

    xa = double(xa);  yb = double(yb);
    if ~isscalar(xa) || ~isscalar(yb)
        [xa, yb] = deal(xa(:), yb(:));
        mode = arrayfun(@(a,b) {paperModeReference(a,b)}, xa, yb);
        mode = [mode{:}].';
        return
    end

    if xa >= 0.28 && yb <= 0.115
        mode = 'FA';
    elseif yb >= max(0.115, 0.205 - 0.40*xa)
        mode = 'AR';
    elseif yb <= 0.145 - 0.50*xa
        mode = 'ST';
    else
        mode = 'CH';
    end
end
