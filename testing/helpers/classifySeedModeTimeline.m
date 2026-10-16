function tl = classifySeedModeTimeline(t, x, opts)
% CLASSIFYSEEDMODETIMELINE  Time-resolved six-part mode label along one trajectory.
%
% classifySeedMode reduces a WHOLE run to one label; this slides a window along
% the trajectory and classifies each window with the same metrics and thresholds,
% so a run whose mode CHANGES (a moving-CoM drop that glides, then spins, then
% dives; a transient before a seed settles) gets a per-sample label for colouring
% a plot or an animation. It replaces the old classifyModeTimeline, which applied
% the retired classifyFlightMode the same way.
%
% The window has to hold a few periods of the motion it is judging: a flip is
% read from the attitude angle's peak-to-peak turns, and a revolution from the
% fitted rate of the span heading. Too short and a slow revolution reads as none.
% A light majority filter over the window labels removes single-window flicker.
%
% INPUTS
%   t    : Nx1 time vector (s).
%   x    : Nx13 state history [r(3) q(4) v(3) omega(3)].
%   opts : struct:
%     .spanLength  seed span (m), REQUIRED: the revolution radius is judged in spans
%     .window      window length (s)                               (default 2.0)
%     .evalDt      spacing of window centres (s)                   (default 0.1)
%     .thresholds  six-part thresholds         (default defaultSeedModeThresholds)
%     .smoothN     majority-filter width in window centres, odd; 1 disables (3)
%     .minSamples  min samples in a window to classify it          (default 8)
%
% OUTPUT (struct tl)
%   .idx       Nx1 mode index per SAMPLE (into .list / .col), nearest-centre mapped
%   .centreT   Cx1 window-centre times (s)
%   .centreIdx Cx1 mode index per window centre (after smoothing)
%   .list      mode names (seedModePalette().names)
%   .col       RGB palette (seedModePalette().colors)

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'spanLength') || isempty(opts.spanLength)
        error('classifySeedModeTimeline:spanLength', ...
              'opts.spanLength is required: the revolution radius is judged in spans.');
    end
    if ~isfield(opts, 'window');     opts.window     = 2.0;  end
    if ~isfield(opts, 'evalDt');     opts.evalDt     = 0.1;  end
    if ~isfield(opts, 'thresholds') || isempty(opts.thresholds)
        opts.thresholds = defaultSeedModeThresholds();
    end
    if ~isfield(opts, 'smoothN');    opts.smoothN    = 3;    end
    if ~isfield(opts, 'minSamples'); opts.minSamples = 8;    end

    pal = seedModePalette();
    t = t(:);   N = numel(t);
    mo = struct('windowStartFrac', 0, 'spanLength', opts.spanLength, ...
                'flipHysteresisTurns', opts.thresholds.flipHysteresisTurns);

    % --- Window centres (centred windows, clipped at the ends) -------------
    tc = (t(1) : opts.evalDt : t(end)).';
    if isempty(tc); tc = t(1); end
    if tc(end) < t(end) - 1e-9; tc(end+1) = t(end); end

    % --- Classify each window ---------------------------------------------
    ci = zeros(numel(tc), 1);              % 0 = not classified yet
    W  = opts.window;
    for j = 1 : numel(tc)
        lo  = max(tc(j) - W/2, t(1));   hi = min(tc(j) + W/2, t(end));
        sel = (t >= lo) & (t <= hi);
        if nnz(sel) < opts.minSamples; continue; end
        try
            lab = classifySeedMode(computeSeedModeMetrics(t(sel), x(sel, :), mo), opts.thresholds);
            ci(j) = find(strcmp(lab, pal.names), 1);
        catch
            % a window too degenerate for the metrics keeps the neighbours' label
        end
    end
    % Unclassified windows take the nearest classified neighbour's label (leading
    % ones the first label); a run with no classifiable window is 'failed'.
    ok = ci > 0;
    if ~any(ok)
        ci(:) = find(strcmp('failed', pal.names), 1);
    elseif ~all(ok)
        ci(~ok) = interp1(find(ok), ci(ok), find(~ok), 'nearest', 'extrap');
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
    tl.list      = pal.names;
    tl.col       = pal.colors;
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
