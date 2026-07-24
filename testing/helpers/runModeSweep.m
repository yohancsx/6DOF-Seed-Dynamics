function results = runModeSweep(sweepSpec, cfg, baselineBsp, baselineSeedParams)
% RUNMODESWEEP  Sweep a nut-position axis, run each drop, classify the mode.
%
% One call = one mode-elicitation experiment: move the nut along an axis over a
% range, integrate each drop, compute trajectory metrics, classify the flight
% mode, print a table (with expected-vs-found), and draw a colour-graded
% overlay of all the trajectories. Reuses runOneSeedCase / computeTrajectoryMetrics
% / classifyFlightMode / plotGroupOverlay.
%
% INPUTS
%   sweepSpec : struct describing the sweep:
%     .name         - label / group name (e.g. 'glide_dive')
%     .axis         - 'chord' | 'span' | 'diag' | 'vertical' (nut offset axis),
%                     or 'C_span' to hold the nut FIXED and sweep the span-force
%                     tuning coefficient instead (see .fixedAxis/.fixedNutFrac).
%     .fracRange    - [lo hi] swept range. For a nut axis: offset as a fraction
%                     of the half-dimension (half-chord for chord/vertical,
%                     half-span for span; diag uses both); values >1 put the nut
%                     OFF the planform. For 'C_span': the C_span values
%                     themselves (0 = span force fully off, since the force
%                     scales linearly with it).
%     .fixedAxis    - ('C_span' sweeps only) nut axis held fixed, e.g. 'diag'
%     .fixedNutFrac - ('C_span' sweeps only) nut offset fraction held fixed
%     .nSweep       - number of increments
%     .expectedMode - char or cellstr of the mode(s) you expect (for the check)
%     .omega0       - (optional) 3x1 initial body angular velocity to seed the
%                     drop (e.g. a small yaw nudge to break symmetry for AR)
%     .q0           - (optional) 4x1 initial orientation quaternion (body->world,
%                     scalar-first). Default identity (level, broadside drop).
%                     Build one with axisAngleToQuat, e.g. a pi/6 tilt about the
%                     body z (spanwise) axis: axisAngleToQuat([0;0;1], pi/6)
%   cfg               : config (geometry, environment, tspan, tolerances,
%                       .metricOpts, .modeThresholds, .plotVisible).
%   baselineBsp       : baseline base-seed-params (nut position gets overridden).
%   baselineSeedParams: pre-built baseline seedParams (unused here since every
%                       case rebuilds, but kept for the runOneSeedCase signature).
%
% OUTPUT
%   results : struct array (one per sweep point) with the runOneSeedCase fields
%             plus .metrics, .mode, .modeInfo.

    fracs = linspace(sweepSpec.fracRange(1), sweepSpec.fracRange(2), sweepSpec.nSweep);
    hc = cfg.chordLength / 2;   % half-chord (body x)
    hs = cfg.spanLength  / 2;   % half-span  (body z)

    if isfield(sweepSpec, 'omega0') && ~isempty(sweepSpec.omega0)
        omega0 = sweepSpec.omega0(:);
    else
        omega0 = [0; 0; 0];
    end

    % Initial orientation: identity (level, broadside) unless the sweep asks for
    % a tilt (e.g. to break the symmetric drop and let tumbling/spiralling start).
    if isfield(sweepSpec, 'q0') && ~isempty(sweepSpec.q0)
        q0 = sweepSpec.q0(:);
        q0 = q0 / norm(q0);           % defensive re-normalisation
    else
        q0 = [1; 0; 0; 0];
    end

    results = [];
    for f = fracs
        % Per-case config copy so any aero override cannot leak between points.
        cfgCase = cfg;

        % --- Nut position (and, for C_span sweeps, the aero override) ------
        if strcmp(sweepSpec.axis, 'C_span')
            % Nut held fixed; the swept quantity is the span-force coefficient.
            pos = nutPosForAxis(sweepSpec.fixedAxis, sweepSpec.fixedNutFrac, hc, hs);
            if isfield(cfgCase, 'aero') && ~isempty(cfgCase.aero)
                cfgCase.aero.C_span = f;      % keep other overrides, replace C_span
            else
                cfgCase.aero = struct('C_span', f);   % partial: rest fall back to defaults
            end
        else
            pos = nutPosForAxis(sweepSpec.axis, f, hc, hs);
        end

        % --- Assemble the case (dropped from rest unless omega0 given) -----
        cs = struct();
        cs.group        = sweepSpec.name;
        cs.label        = sprintf('%s_%+.2f', sweepSpec.name, f);
        cs.bspOverrides = struct('nutPos_t', repmat(pos, 1, numel(baselineBsp.tSamples)));
        cs.postBuildFn  = [];
        cs.q0           = q0;
        cs.v0           = [0; 0; 0];
        cs.omega0       = omega0;
        cs.sweepVal     = f;

        r = runOneSeedCase(cs, baselineSeedParams, baselineBsp, cfgCase);

        % --- Metrics + classification -------------------------------------
        if strcmp(r.status, 'OK')
            r.metrics = computeTrajectoryMetrics(r.t, r.x, cfgCase.metricOpts);
            [r.mode, r.modeInfo] = classifyFlightMode(r.metrics, cfgCase.modeThresholds);
        else
            r.metrics  = [];
            r.mode     = r.status;   % 'FAILED' / 'NONFINITE'
            r.modeInfo = [];
        end

        results = [results; r]; %#ok<AGROW>
    end

    printSweepTable(results, sweepSpec);

    % --- Overlay of all trajectories in the sweep -------------------------
    if ~isfield(cfg, 'plotVisible'); cfg.plotVisible = 'on'; end
    plotGroupOverlay(results, sweepSpec.name, cfg.plotVisible);
end


% =========================================================================
% LOCAL: nut position (body [x;y;z]) for a given offset axis and fraction
% =========================================================================
function pos = nutPosForAxis(axisName, f, hc, hs)
% hc = half-chord (body x reference), hs = half-span (body z reference).
% 'vertical' is the out-of-plane (body y) axis, referenced to the half-chord
% since the plate itself is thin.
    switch axisName
        case 'chord';    pos = [f*hc; 0;    0   ];
        case 'span';     pos = [0;    0;    f*hs];
        case 'diag';     pos = [f*hc; 0;    f*hs];
        case 'vertical'; pos = [0;    f*hc; 0   ];
        otherwise
            error('runModeSweep:badAxis', 'Unknown nut axis "%s".', axisName);
    end
end


% =========================================================================
% LOCAL: print a per-sweep classification table + expected-vs-found summary
% =========================================================================
function printSweepTable(results, sweepSpec)
    expected = sweepSpec.expectedMode;
    if ischar(expected); expected = {expected}; end

    fprintf('\n=== Mode sweep: %s  (axis=%s, expecting: %s) ===\n', ...
            sweepSpec.name, sweepSpec.axis, strjoin(expected, ' | '));
    fprintf('%-8s %-16s %-8s %-8s %-9s %-7s %-6s\n', ...
            'value', 'mode', 'descent', 'glide', 'vSpin', 'cone', 'match');

    nMatch = 0;
    for k = 1:numel(results)
        r = results(k);
        if isempty(r.metrics)
            fprintf('%-8.2f %-16s %-8s %-8s %-9s %-7s %-6s\n', ...
                    r.sweepVal, r.mode, '-', '-', '-', '-', '-');
            continue
        end
        mt   = r.metrics;
        isEx = any(strcmp(r.mode, expected));
        nMatch = nMatch + isEx;
        fprintf('%-8.2f %-16s %-8.2f %-8.2f %-9.2f %-7.1f %-6s\n', ...
                r.sweepVal, r.mode, mt.descentSpeed, mt.glideRatio, ...
                mt.verticalSpinMag, mt.coneAngleDeg, ternary(isEx, 'yes', ''));
    end
    fprintf('--> %d/%d sweep points classified as an expected mode.\n', ...
            nMatch, numel(results));
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end
