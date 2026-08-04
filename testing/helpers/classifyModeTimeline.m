function tl = classifyModeTimeline(t, x, opts)
% CLASSIFYMODETIMELINE  Time-resolved flight-mode label along one trajectory.
%
% classifyFlightMode reduces a WHOLE run to a single mode; this slides a short
% window along the trajectory and classifies each window, so a run whose mode
% CHANGES over time (e.g. a moving-CoM drop that glides, then spins, then dives)
% gets a per-sample mode label suitable for colouring a plot or animation.
%
% Within a short window "converged" is not meaningful (the seed is mid-mode or
% mid-transition), so by default the converged flag is forced true and the
% geometric/rotational rules of classifyFlightMode decide the label. A light
% majority filter over the window-centre labels removes single-window flicker.
%
% INPUTS
%   t    : Nx1 time vector (s).
%   x    : Nx13 state history [r(3) q(4) v(3) omega(3)] (as returned by ode45 on
%          seed6DOFODE; rows = time).
%   opts : (optional) struct:
%     .window         window length in seconds            (default 0.6)
%     .centered       true = window centred on each eval time; false = trailing
%                                                          (default true)
%     .evalDt         spacing of window centres (s)        (default 0.05)
%     .thresholds     classifyFlightMode thresholds        (default defaultModeThresholds)
%     .forceConverged force m.converged=true before classifying (default true)
%     .smoothN        majority-filter width in window-centres, odd; 1 disables
%                                                          (default 3)
%     .minSamples     min samples in a window to classify it (default 8)
%
% OUTPUT (struct tl)
%   .idx      Nx1 mode index per SAMPLE (into .list / .col), nearest-centre mapped
%   .centreT  Cx1 window-centre times (s)
%   .centreIdx Cx1 mode index per window centre (after smoothing)
%   .list     1x10 cell of mode labels (from seedModeColors)
%   .col      10x3 RGB palette (from seedModeColors)

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'window');         opts.window         = 0.6;  end
    if ~isfield(opts, 'centered');       opts.centered       = true; end
    if ~isfield(opts, 'evalDt');         opts.evalDt         = 0.05; end
    if ~isfield(opts, 'thresholds') || isempty(opts.thresholds)
        opts.thresholds = defaultModeThresholds();
    end
    if ~isfield(opts, 'forceConverged'); opts.forceConverged = true; end
    if ~isfield(opts, 'smoothN');        opts.smoothN        = 3;    end
    if ~isfield(opts, 'minSamples');     opts.minSamples     = 8;    end

    t = t(:);   N = numel(t);
    [modeList, modeCol] = seedModeColors();
    undIdx = find(strcmp('undetermined', modeList), 1);

    % --- Window centres ---------------------------------------------------
    tc = (t(1) : opts.evalDt : t(end)).';
    if isempty(tc); tc = t(1); end
    if tc(end) < t(end) - 1e-9; tc(end+1) = t(end); end

    % --- Classify each window ---------------------------------------------
    ci      = zeros(numel(tc), 1);
    lastIdx = undIdx;
    W       = opts.window;
    for j = 1 : numel(tc)
        if opts.centered
            lo = tc(j) - W/2;   hi = tc(j) + W/2;
        else
            lo = tc(j) - W;     hi = tc(j);
        end
        lo = max(lo, t(1));   hi = min(hi, t(end));
        sel = (t >= lo) & (t <= hi);
        if nnz(sel) < opts.minSamples
            ci(j) = lastIdx;   % too short to classify -> hold previous label
            continue
        end
        m = computeTrajectoryMetrics(t(sel), x(sel, :), struct('windowStartFrac', 0));
        if opts.forceConverged; m.converged = true; end
        lab = classifyFlightMode(m, opts.thresholds);
        k = find(strcmp(lab, modeList), 1);
        if isempty(k); k = undIdx; end
        ci(j) = k;   lastIdx = k;
    end

    % --- Light majority filter to de-flicker the centre labels ------------
    if opts.smoothN > 1 && numel(ci) >= opts.smoothN
        ci = majorityFilter(ci, opts.smoothN);
    end

    % --- Map centre labels to every sample (nearest centre) ---------------
    if numel(tc) < 2
        idx = repmat(ci(1), N, 1);
    else
        idx = interp1(tc, ci, t, 'nearest', 'extrap');
    end

    tl.idx       = idx(:);
    tl.centreT   = tc;
    tl.centreIdx = ci;
    tl.list      = modeList;
    tl.col       = modeCol;
end


% =========================================================================
% LOCAL: centred majority (statistical-mode) filter over integer labels
% =========================================================================
function y = majorityFilter(v, w)
    n = numel(v);   y = v;   h = floor(w / 2);
    for i = 1 : n
        lo = max(1, i - h);   hi = min(n, i + h);
        y(i) = mode(v(lo:hi));   % most-frequent label in the neighbourhood
    end
end
