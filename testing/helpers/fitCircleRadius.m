function [R, valid, centre] = fitCircleRadius(u, v)
% FITCIRCLERADIUS  Algebraic (Kasa) circle fit to a 2D track; radius + centre.
%
% Used for the horizontal ground track of a descending seed: a genuine helix
% fits a finite radius, a near-straight track does not. Promoted out of
% computeTrajectoryMetrics (where it was a local function) unchanged, so that
% computePaperMetrics can also reach the CENTRE -- which is what the orbital
% revolution rate is measured about.
%
% INPUTS
%   u, v : the two horizontal coordinates of the track (world X and Z here).
% OUTPUTS
%   R      : fitted radius (m), or Inf when the fit is not meaningful.
%   valid  : true when the points really do lie near a circle (see below).
%   centre : 2x1 fitted centre [u_c; v_c], or [NaN; NaN] when ~valid.

    centre = [NaN; NaN];
    u = u(:);  v = v(:);   n = numel(u);
    if n < 5
        R = Inf; valid = false; return
    end

    % Spatial extent of the horizontal track. A genuine helix curves back on
    % itself, so its fitted radius is comparable to this extent; a near-straight
    % track has ~zero extent-normalised curvature and the fit below is
    % ill-conditioned, throwing R off to ~1e9.
    extent = hypot(max(u) - min(u), max(v) - min(v));
    if extent < eps
        R = Inf; valid = false; return
    end

    A = [u, v, ones(n, 1)];
    b = u.^2 + v.^2;
    % A near-collinear (straight) track makes A rank-deficient; solve quietly and
    % judge the result by the checks below rather than trusting the warning.
    ws  = warning('off', 'MATLAB:rankDeficientMatrix');
    sol = A \ b;
    warning(ws);

    uc = sol(1) / 2;
    vc = sol(2) / 2;
    R  = sqrt(max(sol(3) + uc^2 + vc^2, 0));

    % Valid only if the points actually lie near a circle: small radial spread
    % relative to R, AND a radius that is not absurdly large versus the track
    % extent (a straight path fits R >> extent). Otherwise it is not a helix:
    % report R = Inf so downstream code never sees a spurious finite radius.
    resid = std(hypot(u - uc, v - vc));
    valid = isfinite(R) && R > 0 && resid < 0.5 * R && R < 50 * extent;
    if valid
        centre = [uc; vc];
    else
        R = Inf;
    end
end
