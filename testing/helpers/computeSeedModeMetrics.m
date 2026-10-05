function sm = computeSeedModeMetrics(t, x, opts)
% COMPUTESEEDMODEMETRICS  Metrics for classifySeedMode (the six-part classifier).
%
% Everything computePaperMetrics returns, plus the quantities the six-part scheme
% in testing/classifier/MODE_DEFINITIONS.md is defined on. Same settled window.
%
% NEW QUANTITIES
%   psi (flip angle), measured from ATTITUDE, not integrated from rates.
%       The plate normal's angle about the span axis, referenced to "up":
%           e1 = unit( Yhat - (Yhat.shat) shat )     up, perpendicular to the span
%           e2 = shat x e1
%           psi = atan2( nhat.e2, nhat.e1 ),  unwrapped
%       A steady revolution leaves psi constant (whatever its rate or radius); a
%       tumble advances it 2*pi per turn; a flutter makes it oscillate. Nothing
%       here needs the revolution rate, so a mis-fitted or absent helix cannot
%       leak into the flip measurement. Undefined only when the span is
%       vertical -- the end-on case, which the classifier catches first.
%   flipPPTurns       (max psi - min psi) / 2pi                    turns
%   flipMaxTurnsOneWay largest one-way leg of psi between reversals turns
%   flipReversals     direction reversals of psi, with a hysteresis of
%                     opts.flipHysteresisTurns (default 0.25 turn): psi must
%                     retreat that far from its running extreme to count
%   flipReversalCV    std/mean of the time between reversals (NaN if < 3
%                     reversals) -- small = periodic (segmented tumbling),
%                     large = irregular
%   spanHeadingRate   slope of the span axis's horizontal heading, measured
%                     about +Y (same sense as omega.Yhat)            rad/s
%                     -- the REVOLUTION rate, from attitude. A tumble or
%                     flutter about the span leaves the heading fixed
%                     however the span is tilted; the body spin about the
%                     vertical does not (a tilted tumble projects onto it).
%                     Undefined for a vertical span (end-on, caught first).
%   absSpanTiltDeg    mean |theta|, span axis from horizontal       deg
%   radiusPerSpan     helixRadius / opts.spanLength  (NaN if no valid helix)
%   spanAlongVel      mean |shat . vhat|: how much the seed moves along its span
%
% INPUTS
%   t, x : ode45 outputs (Nx1, Nx13).
%   opts : .windowStartFrac, .convergeTol (as computeTrajectoryMetrics), and
%          .spanLength (m, REQUIRED -- the radius is judged in spans),
%          .flipHysteresisTurns (optional, default 0.25).
%
% OUTPUT
%   sm : the computePaperMetrics struct with the fields above added.

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'spanLength') || isempty(opts.spanLength) || ~(opts.spanLength > 0)
        error('computeSeedModeMetrics:noSpan', ...
              'opts.spanLength is required: the revolution radius is judged in spans.');
    end
    if ~isfield(opts, 'windowStartFrac');     opts.windowStartFrac     = 0.5;  end
    if ~isfield(opts, 'flipHysteresisTurns'); opts.flipHysteresisTurns = 0.25; end

    sm = computePaperMetrics(t, x, opts);

    t = t(:);
    win = t >= opts.windowStartFrac * t(end);
    if nnz(win) < 10; win = true(numel(t), 1); end
    iw = find(win);   nw = numel(iw);   tw = t(iw);

    % --- psi from attitude, span tilt, motion along the span ----------------
    psi = zeros(nw,1);   absTh = zeros(nw,1);   along = zeros(nw,1);   az = zeros(nw,1);
    Y = [0;1;0];
    for kk = 1:nw
        k = iw(kk);
        q = x(k,4:7).';   R = quatToRotm(q / norm(q));
        s = R*[0;0;1];   n = R*[0;1;0];
        e1 = Y - (Y.'*s)*s;   ne = norm(e1);
        if ne > 1e-9
            e1 = e1 / ne;   e2 = cross(s, e1);
            psi(kk) = atan2(n.'*e2, n.'*e1);
        elseif kk > 1
            psi(kk) = psi(kk-1);          % span exactly vertical: hold
        end
        absTh(kk) = abs(asind(max(-1, min(1, s(2)))));
        az(kk)    = atan2(-s(3), s(1));   % span heading, measured about +Y
        v = x(k,8:10).';   sp = norm(v);
        if sp > 0;  along(kk) = abs(s.'*v) / sp;  end
    end
    psi = unwrap(psi);

    % --- flip legs and reversals (zig-zag with hysteresis, in angle) ---------
    h = 2*pi*opts.flipHysteresisTurns;
    piv = psi(1);  hi = psi(1);  lo = psi(1);  dir = 0;  legs = [];  revT = [];
    for kk = 2:nw
        p = psi(kk);   hi = max(hi, p);   lo = min(lo, p);
        switch dir
            case 0
                if hi - piv > h;      dir = 1;   lo = p;
                elseif piv - lo > h;  dir = -1;  hi = p;
                end
            case 1                            % rising: reversal when p falls h below hi
                if hi - p > h
                    legs(end+1) = hi - piv;   piv = hi;   revT(end+1) = tw(kk); %#ok<AGROW>
                    dir = -1;   lo = p;
                end
            case -1                           % falling: reversal when p rises h above lo
                if p - lo > h
                    legs(end+1) = piv - lo;   piv = lo;   revT(end+1) = tw(kk); %#ok<AGROW>
                    dir = 1;    hi = p;
                end
        end
    end
    if dir == 1;      legs(end+1) = hi - piv;
    elseif dir == -1; legs(end+1) = piv - lo;
    else;             legs(end+1) = max(hi - piv, piv - lo);
    end

    sm.flipPPTurns        = (max(psi) - min(psi)) / (2*pi);
    sm.flipMaxTurnsOneWay = max(legs) / (2*pi);
    sm.flipReversals      = numel(revT);
    if numel(revT) >= 3
        dT = diff(revT);
        sm.flipReversalCV = std(dT) / mean(dT);
    else
        sm.flipReversalCV = NaN;
    end
    pAz = polyfit(tw - tw(1), unwrap(az), 1);
    sm.spanHeadingRate = pAz(1);
    sm.absSpanTiltDeg = mean(absTh);
    sm.spanAlongVel   = mean(along);
    if sm.helixValid
        sm.radiusPerSpan = sm.helixRadius / opts.spanLength;
    else
        sm.radiusPerSpan = NaN;
    end
end
