function pm = computePaperMetrics(t, x, opts)
% COMPUTEPAPERMETRICS  Trajectory descriptors in the vocabulary of Hou et al. (2025).
%
% Everything computeTrajectoryMetrics returns, PLUS the four quantities that
% Hou, Zhang, Li, Jia & Huang, "Aerodynamic significance of mass distribution on
% diverse samara descent behaviors", Commun. Eng. 4, 129 (2025) actually report --
% so their published numbers can be compared to ours directly, without
% reinterpreting either side.
%
% WHY THESE ARE NEW (and not just renames of what we already had):
%
%   theta (spanAxisTiltDeg)  Their CONE ANGLE is the angle between the plate's
%       LONGER axis and the horizontal plane. Our coneAngleDeg is the angle of the
%       plate NORMAL from vertical. Those are different quantities: a plate can
%       hold its long axis at a steady tilt while its normal swings through 180
%       degrees (that is exactly what spiral tumbling does). Their values:
%       AR -11.4 +/- 4.2 deg, continuous ST -38.2 +/- 2.3, segmented ST -38.3 +/- 13.0.
%       SIGN CONVENTION HERE: +z is the nut side of the span, so theta < 0 means
%       the nut-side tip is BELOW horizontal -- which matches the sign they report.
%
%   phi' (revolutionRate)  Their revolution rate is ORBITAL: how fast the plate
%       travels around the vertical axis of its helix. Our verticalSpinMag is the
%       body angular velocity projected on the world vertical, which is a different
%       thing entirely for a tumbling plate. This is measured from the horizontal
%       ground track about the fitted circle centre, so it is the orbital rate by
%       construction. Signed: + is counter-clockwise seen from above.
%
%   psi' (selfRotationRate)  Their SELF-ROTATION about the long axis -- the
%       tumbling proper. This is NOT the body-frame omega_z, and the difference
%       matters: a plate that merely REVOLVES at a tilt already has a constant
%       nonzero omega_z, because the revolution vector projects onto the tilted
%       span axis. Decomposing omega = phi'*Yhat + psi'*shat and projecting on
%       shat gives
%             psi' = omega_z - phi' * sin(theta),      sin(theta) = shat . Yhat
%       so the revolution's projection has to be removed before anything is called
%       tumbling. Their Fig. 3f plots <|phi'|> and <|psi'|> as separate
%       quantities, which is exactly this decomposition.
%       (Using omega_z raw made a pure 11-degree-tilt revolution at 40 rad/s read
%       as 7.9 rad/s of "continuous tumbling" -- it would have been labelled CST
%       instead of AR. testPaperMetrics case A is that exact trajectory.)
%
%   tumblesPerRev  Their "roughly seven tumbling rotations per spiral revolution"
%       for the continuous ST mode. = <|psi'|> / |phi'|.
%
%   tumbleReversals  What makes SEGMENTED spiral tumbling segmented: at each
%       turning point the tumbling direction reverses (their Fig. 3e, CCW -> CW).
%       Counted with a Schmitt trigger on a smoothed psi', so noise near zero
%       cannot manufacture reversals.
%
% Optional normalisations, computed when the corresponding opts field is given:
%   helixRadiusPerL = helixRadius / opts.refLength   (they quote R ~ 1.25 L and 1.89 L)
%   VdScaled        = descentSpeed / sqrt(sigma*g/rho)   (their scaling law, eq. 2;
%                     the prefactor is O(1), so this number says how far off our
%                     descent speed is in a geometry-independent way)
%
% INPUTS
%   t    : Nx1 time vector (s).
%   x    : Nx13 state history [r(3) q(4) v(3) omega(3)], world r/v, body omega.
%   opts : (optional) struct. Passed through to computeTrajectoryMetrics
%          (.windowStartFrac, .convergeTol), plus:
%          .refLength  - length to normalise the helix radius by (m), e.g. the span.
%          .sigma      - average surface density, total mass / planform area (kg/m^2).
%          .g, .rhoFluid - to scale with (sigma*g/rho).
%
% OUTPUT
%   pm : the computeTrajectoryMetrics struct with these fields added:
%        .spanAxisTiltDeg .spanAxisTiltStd .revolutionRate .revolutionRateMag
%        .selfRotationRate .selfTumbleFrac .maxTurnsOneWay .tumbleRate
%        .tumblesPerRev .tumbleReversals .helixCentre .helixRadiusPerL .VdScaled
%        (.tumbleRate is the raw mean |omega_z|, kept only for reference; the
%         tumbling quantity is .selfRotationRate.)

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'windowStartFrac'); opts.windowStartFrac = 0.5;  end
    if ~isfield(opts, 'convergeTol');     opts.convergeTol     = 0.20; end

    pm = computeTrajectoryMetrics(t, x, opts);

    t = t(:);
    N = numel(t);
    r = x(:, 1:3);   q = x(:, 4:7);   w = x(:, 11:13);

    % Same steady-state window as computeTrajectoryMetrics, by the same rule.
    win = t >= opts.windowStartFrac * t(end);
    if nnz(win) < 10; win = true(N, 1); end
    iw = find(win);
    nw = numel(iw);

    % --- theta: the LONG axis (body z) from the horizontal plane -----------
    tilt    = zeros(nw, 1);      % deg
    tiltSin = zeros(nw, 1);      % shat . Yhat, used to remove the revolution below
    for kk = 1:nw
        R = quatToRotm(q(iw(kk), :).');
        sWorld  = R * [0; 0; 1];                       % span axis in world, +z = nut side
        tiltSin(kk) = max(-1, min(1, sWorld(2)));
        tilt(kk)    = asind(tiltSin(kk));              % + = nut-side tip above horizontal
    end
    pm.spanAxisTiltDeg = mean(tilt);
    pm.spanAxisTiltStd = std(tilt);

    % --- phi': orbital rate about the fitted helix axis --------------------
    X = r(iw, 1);   Z = r(iw, 3);
    [~, hValid, centre] = fitCircleRadius(X, Z);
    pm.helixCentre = centre;
    if hValid
        % Angle measured about +Y (right-handed, world Y up): a rotation by +phi
        % about +Y carries the X axis toward -Z, so the angle is atan2(-dZ, dX).
        % This is the SAME sense as omega . Yhat, which the psi' decomposition
        % below requires. (Fixed 2026-10: this used atan2(+dZ, dX), the angle
        % about -Y, so phi' had the opposite sign to the body rotation and
        % psi' = omega_z - phi'*sin(theta) DOUBLED a steady revolution instead of
        % cancelling it -- every steady revolver read as continuous tumbling.
        % testPaperMetrics' synthetic orbit had the matching flip, so it passed.)
        phi = unwrap(atan2(-(Z - centre(2)), X - centre(1)));
        p   = polyfit(t(iw), phi, 1);
        pm.revolutionRate = p(1);
    else
        pm.revolutionRate = NaN;       % not a helix: no orbital rate to measure
    end
    pm.revolutionRateMag = abs(pm.revolutionRate);

    % --- psi': self-rotation about the long axis, and its reversals --------
    % omega_z alone is NOT the tumbling rate: subtract the revolution's
    % projection on the tilted span axis first (see header).
    wz = w(iw, 3);
    pm.tumbleRate = mean(abs(wz));            % raw |omega_z|, kept for reference
    if isfinite(pm.revolutionRate)
        psiDot = wz - pm.revolutionRate * tiltSin;
    else
        psiDot = wz;                          % not revolving: nothing to remove
    end
    pm.selfRotationRate = mean(abs(psiDot));
    pm.selfTumbleFrac   = abs(trapz(t(iw), psiDot)) / max(trapz(t(iw), abs(psiDot)), eps);
    if isfinite(pm.revolutionRateMag) && pm.revolutionRateMag > 0
        pm.tumblesPerRev = pm.selfRotationRate / pm.revolutionRateMag;
    else
        pm.tumblesPerRev = NaN;
    end
    [pm.tumbleReversals, pm.maxTurnsOneWay] = analysePsi(psiDot, t(iw), nw);

    % --- optional normalisations -------------------------------------------
    if isfield(opts, 'refLength') && ~isempty(opts.refLength) && opts.refLength > 0
        pm.helixRadiusPerL = pm.helixRadius / opts.refLength;
    else
        pm.helixRadiusPerL = NaN;
    end
    if all(isfield(opts, {'sigma', 'g', 'rhoFluid'}))
        pm.VdScaled = pm.descentSpeed / sqrt(opts.sigma * opts.g / opts.rhoFluid);
    else
        pm.VdScaled = NaN;
    end
end


% =========================================================================
% LOCAL: reversals of the tumbling direction, and the largest one-way rotation
% =========================================================================
function [n, maxTurns] = analysePsi(psiDot, tw, nw)
% Smooth psi' over ~2% of the window, then run a Schmitt trigger at 25% of the
% peak, so a signal that merely grazes zero (or noise about zero) is not read as a
% reversal. Returns:
%   n         reversals of the tumbling direction. Continuous ST -> 0; segmented
%             -> one per turning point (their Fig. 3e).
%   maxTurns  the largest NET rotation accumulated in a single direction, in
%             turns. This is what separates TUMBLING from FLUTTERING, and it has
%             to be measured per-segment rather than over the whole window:
%             segmented ST reverses, so its net rotation over the window is near
%             zero even though it is unmistakably tumbling between turning
%             points. A tumbler completes whole revolutions one way (>= 1 turn);
%             a flutterer oscillates through a limited angle and never does.
    k  = max(3, round(0.02 * nw));
    ws = movmean(psiDot, k);
    pk = max(abs(ws));
    if pk <= 0;  n = 0;  maxTurns = 0;  return;  end
    thr = 0.25 * pk;

    n = 0;  state = 0;  segStart = 1;  maxTurns = 0;
    for i = 1:numel(ws)
        if     ws(i) >  thr;  s =  1;
        elseif ws(i) < -thr;  s = -1;
        else                  s =  0;      % inside the deadband: hold
        end
        if s ~= 0
            if state ~= 0 && s ~= state
                n = n + 1;
                maxTurns = max(maxTurns, segTurns(psiDot, tw, segStart, i));
                segStart = i;
            end
            state = s;
        end
    end
    maxTurns = max(maxTurns, segTurns(psiDot, tw, segStart, numel(ws)));
end

function turns = segTurns(psiDot, tw, i1, i2)
    if i2 - i1 < 2;  turns = 0;  return;  end
    turns = abs(trapz(tw(i1:i2), psiDot(i1:i2))) / (2*pi);
end
