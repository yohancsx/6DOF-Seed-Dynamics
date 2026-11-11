function out = runRealSeedSet3D(cfg, outDir)
% RUNREALSEEDSET3D  Fly the measured seeds and compare them with their tracked drops.
%
%   out = runRealSeedSet3D()                    all seeds, defaults, output to P.misc
%   out = runRealSeedSet3D(cfg, outDir)         as stage G of runSeedTestSuite3D
%
% VALIDATION ONLY. Nothing is fitted to the drops: the seeds are the CT-fitted
% plates and mass models (loadRealSeeds), flown with the model as it stands.
%
% WHAT IT DOES, per seed and release
%   1. Builds the seed (buildSeedParams on loadRealSeeds' stripGeom). For real seeds
%      the per-strip added-mass form is the default here (enableAddedMass3D = true;
%      the tips carry up to ~26 deg dihedral), unless cfg sets it.
%   2. Releases it from rest at cfg.realDropHeight (10 m) in each attitude of
%      cfg.realReleasesDeg about the SPAN axis (0 = flat, 135), no spin, and
%      integrates to the ground (ode45 event), capped at cfg.realTimeCap.
%   3. Labels the motion with the six-part classifier (the simulation has attitude).
%   4. Measures the simulated CoM path, resampled at the camera rate and put in the
%      lab frame (z up), and the tracked drop with ONE function, computePathMetrics,
%      over the SAME band of height fallen: from cfg.realWindowStart to where the
%      shorter of the two ends (some seeds glide out of the camera view).
% The real release attitude is unknown, so both releases are reported side by side
% and neither is chosen as "the" match; their difference is the release sensitivity.
%
% LAB FRAME. Simulation world: Y up. Lab: z up. lab = [X, -Z, Y] (right-handed, so
% the sense of a helix -- its chirality -- is preserved).
%
% CFG FIELDS (defaults in brackets)
%   realSeedFile [P.seedMassModel]   trackFile [P.seedTracks]   trackFps [60]
%   realSeeds [all]   realCalibration ['both']   realSubStrips [1]
%   realReleasesDeg [0 135]   realDropHeight [10]   realTimeCap [60]
%   realWindowStart [1.5]   realSeedVideos [false]   realPlots [true]
%   plus the usual rhoFluid, g, odeRelTol, odeAbsTol, metricOpts,
%   seedModeThresholds, and any physics switch overrides (passed to buildSeedParams).
%
% OUTPUTS
%   files in outDir: real_seeds_comparison.csv (one row per seed x release),
%                    real_seeds_3D.mat, summary_real_seeds.txt, and figures
%                    real_seeds_descent.png, real_seeds_paths_side.png,
%                    real_seeds_paths_top.png, real_seeds_labels.png.
%   out : .T (the comparison table), .agg (headline numbers), .lines (summary text),
%         .tracked (the numbers for the suite's tracking history), .runs.

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  P = seedPaths(here);
    if nargin < 1 || isempty(cfg); cfg = struct(); end
    if nargin < 2 || isempty(outDir)
        outDir = fullfile(P.misc, 'Real Seeds', datestr(now, 'yyyy-mm-dd_HHMMSS')); %#ok<TNOW1,DATST>
    end
    d = struct('realSeedFile', P.seedMassModel, 'trackFile', P.seedTracks, 'trackFps', 60, ...
               'realSeeds', [], 'realCalibration', 'both', 'realSubStrips', 1, ...
               'realReleasesDeg', [0 135], 'realDropHeight', 10, 'realTimeCap', 60, ...
               'realWindowStart', 1.5, 'realSeedVideos', false, 'realPlots', true, ...
               'rhoFluid', 1.225, 'g', 9.81, 'odeRelTol', 1e-6, 'odeAbsTol', 1e-8, ...
               'metricOpts', struct('windowStartFrac', 0.5, 'convergeTol', 0.20), ...
               'seedModeThresholds', defaultSeedModeThresholds());
    for f = fieldnames(d).'
        if ~isfield(cfg, f{1}) || (isempty(cfg.(f{1})) && ~strcmp(f{1}, 'realSeeds'))
            cfg.(f{1}) = d.(f{1});
        end
    end
    if ~isfield(cfg, 'enableAddedMass3D'); cfg.enableAddedMass3D = true; end
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    seeds  = loadRealSeeds(cfg.realSeedFile, struct('seeds', cfg.realSeeds, ...
                 'calibration', cfg.realCalibration, 'subStrips', cfg.realSubStrips));
    tracks = loadTrackedSeeds(cfg.trackFile, struct('fps', cfg.trackFps));
    tNums  = [tracks.seedNum];
    noModel = [];                                 % tracked seeds with no mass model
    if isempty(cfg.realSeeds); noModel = setdiff(tNums, [seeds.seedNum]); end

    rel  = cfg.realReleasesDeg(:).';
    jobs = [repelem((1:numel(seeds)).', numel(rel)), repmat((1:numel(rel)).', numel(seeds), 1)];
    nJ   = size(jobs, 1);
    runs = cell(nJ, 1);
    th   = cfg.seedModeThresholds;
    tic;
    parfor jj = 1:nJ
        runs{jj} = flyOne(seeds(jobs(jj,1)), rel(jobs(jj,2)), tracks, tNums, cfg, th);
    end
    fprintf('  real seeds: %d runs in %.0f s\n', nJ, toc);
    runs = [runs{:}];

    % --- comparison table --------------------------------------------------------
    T = struct2table(arrayfun(@rowOf, runs));
    writetable(T, fullfile(outDir, 'real_seeds_comparison.csv'));
    agg = aggregate(T, rel);
    lines = summaryLines(T, agg, rel, cfg, noModel);
    fid = fopen(fullfile(outDir, 'summary_real_seeds.txt'), 'w');  fprintf(fid, '%s\n', lines);  fclose(fid);
    runsLite = rmfield(runs, {'x'});                 % the full state history is large
    save(fullfile(outDir, 'real_seeds_3D.mat'), 'T', 'agg', 'runsLite', 'cfg', '-v7.3');
    if cfg.realPlots
        try
            plotComparison(T, runs, tracks, tNums, rel, outDir);
        catch err
            warning('runRealSeedSet3D:plot', 'plotting failed: %s', err.message);
        end
    end
    if cfg.realSeedVideos
        vdir = fullfile(outDir, 'videos');  if ~exist(vdir, 'dir'); mkdir(vdir); end
        for r = runs
            if isempty(r.x); continue; end
            sp = buildSeedParams(seeds([seeds.seedNum] == r.seedNum).bsp, cfg);
            animateCaseVideo(r.t, r.x, sp, sprintf('seed%d_rel%03d', r.seedNum, r.releaseDeg), vdir, ...
                             struct('modeLabel', char(r.simLabel)));
        end
    end
    out = struct('T', T, 'agg', agg, 'lines', lines, 'runs', runsLite, 'tracked', agg.tracked);
end

% =========================================================================
function r = flyOne(seed, relDeg, tracks, tNums, cfg, th)
    r = struct('seedNum', seed.seedNum, 'releaseDeg', relDeg, 'meta', seed.meta, 'ok', false, ...
               'err', "", 'reachedGround', false, 't', [], 'x', [], 'tq', [], 'pos', [], ...
               'simLabel', "failed", 'simParts', "", 'physical', false, 'sm', [], ...
               'pmSim', [], 'pmExp', [], 'hasTrack', false, 'trackPartial', false, 'window', [NaN NaN]);
    try
        sp  = buildSeedParams(seed.bsp, cfg);
        rhs = seedRHS(sp);
        if relDeg == 0; q0 = [1;0;0;0]; else; q0 = axisAngleToQuat([0;0;1], deg2rad(relDeg)); end
        x0  = [0; cfg.realDropHeight; 0; q0; zeros(6,1)];
        opt = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol, 'Events', @ground);
        [t, x, te] = ode45(@(tt,xx) rhs(tt,xx,sp), [0 cfg.realTimeCap], x0, opt);
        r.t = t;  r.x = x;  r.reachedGround = ~isempty(te);
        so = cfg.metricOpts;  so.spanLength = seed.meta.span;  so.flipHysteresisTurns = th.flipHysteresisTurns;
        sm = computeSeedModeMetrics(t, x, so);
        [lbl, p] = classifySeedMode(sm, th);
        r.sm = sm;  r.simLabel = string(lbl);  r.physical = logical(p.physical);
        r.simParts = string(sprintf('%s/%s/%s/%s/%s', p.flip, p.revolution, p.attitude, p.path, p.regularity));
        % the CoM path at the camera rate, in the lab frame (z up)
        tq = (0:1/cfg.trackFps:t(end)).';
        X  = interp1(t, x(:,1:3), tq);
        pos = [X(:,1), -X(:,3), X(:,2)];
        r.tq = tq;  r.pos = pos;
        simFall = pos(1,3) - min(pos(:,3));
        k = find(tNums == seed.seedNum, 1);
        if ~isempty(k)
            tr = tracks(k);  r.hasTrack = true;  r.trackPartial = tr.partial;
            win = [cfg.realWindowStart, min(simFall, tr.fallen)];
            r.window = win;
            r.pmExp = computePathMetrics(tr.t, tr.pos, struct('window', win));
        else
            win = [cfg.realWindowStart, simFall];  r.window = win;
        end
        r.pmSim = computePathMetrics(tq, pos, struct('window', win));
        r.ok = true;
    catch err
        r.err = string(err.message);
    end
end

function [value, isterminal, direction] = ground(~, x)
    value = x(2);  isterminal = 1;  direction = -1;      % world Y reaches 0
end

function row = rowOf(r)
    m = r.meta;
    row = struct('seed', r.seedNum, 'release_deg', r.releaseDeg, 'calibration', m.calibration, ...
        'ok', r.ok, 'reachedGround', r.reachedGround, 'simLabel', r.simLabel, 'simParts', r.simParts, ...
        'physical', r.physical, 'hasTrack', r.hasTrack, 'trackPartial', r.trackPartial, ...
        'window0_m', r.window(1), 'window1_m', r.window(2));
    for side = {'Sim', 'Exp'}
        pm = r.(['pm' side{1}]);
        if isempty(pm) || ~pm.valid
            v = nan(1,8);  cls = "none";
        else
            h = pm.helix;
            v = [pm.descentSpeed, pm.horizSpeed, pm.driftSpeed, pm.glideDeg, h.radius, h.period, h.turns, h.valid];
            cls = pm.pathClass;
        end
        s = lower(side{1});
        row.(['Vd_' s]) = v(1);       row.(['Vh_' s]) = v(2);   row.(['drift_' s]) = v(3);
        row.(['glide_' s]) = v(4);    row.(['helixR_' s]) = v(5);
        row.(['helixT_' s]) = v(6);   row.(['helixTurns_' s]) = v(7);
        row.(['helixValid_' s]) = v(8);   row.(['class_' s]) = cls;
    end
    row.VdRatio = row.Vd_sim / row.Vd_exp;
    row.glideErr_deg = row.glide_sim - row.glide_exp;
    row.classMatch = row.class_sim == row.class_exp;
    row.massKg = m.M_total;  row.nutFrac = m.nutFrac;  row.span_m = m.span;  row.AR = m.AR;
    row.wingLoading = m.wingLoading;  row.twistRange_deg = m.twistRangeDeg;
    row.confidence = m.confidence;  row.err = r.err;
end

function agg = aggregate(T, rel)
    ok = T.ok & T.hasTrack & isfinite(T.Vd_sim) & isfinite(T.Vd_exp);
    agg.nRuns = height(T);  agg.nFailed = nnz(~T.ok);  agg.nCompared = nnz(ok);
    agg.nonPhysical = mean(~T.physical(T.ok));
    agg.notGround = nnz(T.ok & ~T.reachedGround);
    agg.byRelease = struct('release', {}, 'VdLogErr', {}, 'VdBias', {}, 'VdRank', {}, ...
                           'glideErr', {}, 'classAgree', {}, 'helixAgree', {}, 'nonPhysical', {});
    for a = rel
        w = ok & T.release_deg == a;
        agg.byRelease(end+1) = struct('release', a, ...
            'VdLogErr', median(abs(log(T.VdRatio(w)))), 'VdBias', median(T.VdRatio(w)), ...
            'VdRank', rankCorr(T.Vd_sim(w), T.Vd_exp(w)), ...
            'glideErr', median(abs(T.glideErr_deg(w))), 'classAgree', mean(T.classMatch(w)), ...
            'helixAgree', mean(T.helixValid_sim(w) == T.helixValid_exp(w)), ...
            'nonPhysical', mean(~T.physical(T.ok & T.release_deg == a))); %#ok<AGROW>
    end
    % release sensitivity: same seed, different release -> different label / speed
    seedsAll = unique(T.seed(T.ok));  chg = false(size(seedsAll));  dv = nan(size(seedsAll));
    for i = 1:numel(seedsAll)
        w = T.ok & T.seed == seedsAll(i);
        chg(i) = numel(unique(T.simLabel(w))) > 1;
        v = T.Vd_sim(w);  if numel(v) > 1; dv(i) = (max(v) - min(v)) / mean(v); end
    end
    agg.releaseLabelChange = mean(chg);
    agg.releaseVdSpread = median(dv, 'omitnan');
    b = agg.byRelease;
    agg.tracked = struct('real_VdLogErr', mean([b.VdLogErr]), 'real_VdBias', mean([b.VdBias]), ...
        'real_VdRank', mean([b.VdRank]), 'real_glideErr', mean([b.glideErr]), ...
        'real_classAgree', mean([b.classAgree]), 'real_nonPhysical', agg.nonPhysical, ...
        'real_releaseLabelChange', agg.releaseLabelChange, 'real_nCompared', agg.nCompared);
end

function rho = rankCorr(a, b)
% Spearman rank correlation (no toolbox): Pearson correlation of the ranks.
    a = a(:);  b = b(:);  n = numel(a);
    if n < 3; rho = NaN; return; end
    rk = @(v) localRank(v);
    ra = rk(a);  rb = rk(b);
    c  = corrcoef(ra, rb);  rho = c(1,2);
end
function r = localRank(v)
    [~, o] = sort(v);  r = zeros(size(v));  r(o) = 1:numel(v);
    [u, ~, g] = unique(v);                         % average ties
    if numel(u) < numel(v); r = accumarray(g, r, [], @mean); r = r(g); end
end

function L = summaryLines(T, agg, rel, cfg, noModel)
    L = strings(0,1);
    L(end+1) = sprintf('REAL SEEDS vs TRACKED DROPS  (calibration %s, %d runs: %d seeds x releases %s deg, from %g m)', ...
        cfg.realCalibration, agg.nRuns, numel(unique(T.seed)), mat2str(rel), cfg.realDropHeight);
    L(end+1) = sprintf('  failed %d, did not reach the ground %d, compared with a track %d;  non-physical labels %.0f%%', ...
        agg.nFailed, agg.notGround, agg.nCompared, 100*agg.nonPhysical);
    for b = agg.byRelease
        L(end+1) = sprintf(['  release %3d deg: descent sim/real median %.2f (|log| %.2f), rank corr %+.2f;  ' ...
            'glide err %.1f deg;  path class agree %.0f%%, helix agree %.0f%%;  non-physical %.0f%%'], ...
            b.release, b.VdBias, b.VdLogErr, b.VdRank, b.glideErr, 100*b.classAgree, 100*b.helixAgree, ...
            100*b.nonPhysical); %#ok<AGROW>
    end
    L(end+1) = sprintf('  release sensitivity: label changes with release for %.0f%% of seeds; descent spread median %.0f%%', ...
        100*agg.releaseLabelChange, 100*agg.releaseVdSpread);
    if ~isempty(noModel); L(end+1) = sprintf('  tracked but no mass model (not flown): seed %s', mat2str(noModel)); end
    L(end+1) = "  seed rel | sim label            sim class   real class  | Vd sim  real  | glide sim real | helix R sim real  T sim real";
    for i = 1:height(T)
        r = T(i,:);
        flag = "";  if r.trackPartial; flag = " (track partial)"; end
        if ~r.ok; flag = " FAILED: " + r.err; end
        L(end+1) = sprintf('  %4d %3d | %-20s %-11s %-11s | %5.2f %5.2f | %5.1f %5.1f | %5.2f %5.2f  %5.2f %5.2f%s', ...
            r.seed, r.release_deg, r.simLabel, r.class_sim, r.class_exp, r.Vd_sim, r.Vd_exp, r.glide_sim, ...
            r.glide_exp, r.helixR_sim, r.helixR_exp, r.helixT_sim, r.helixT_exp, flag); %#ok<AGROW>
    end
end

function plotComparison(T, runs, tracks, tNums, rel, outDir)
    col = lines(max(numel(rel), 2));
    % --- descent speed and glide: simulated vs measured -----------------------
    f = figure('Visible', 'off', 'Position', [100 100 1200 420]);
    tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
    names = {'Vd', 'glide', 'helixR'};  labs = {'descent speed (m/s)', 'glide angle (deg)', 'helix radius (m)'};
    for p = 1:3
        ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
        a = T.([names{p} '_exp']);  b = T.([names{p} '_sim']);
        lim = [0, max([a; b; 0.1], [], 'omitnan')*1.1];
        plot(ax, lim, lim, 'k:');
        for i = 1:numel(rel)
            w = T.release_deg == rel(i);
            scatter(ax, a(w), b(w), 28, col(i,:), 'filled', 'DisplayName', sprintf('release %d deg', rel(i)));
        end
        xlim(ax, lim);  ylim(ax, lim);  axis(ax, 'square');
        xlabel(ax, ['measured ' labs{p}]);  ylabel(ax, ['simulated ' labs{p}]);
        if p == 1; legend(ax, 'Location', 'northwest'); end
    end
    title(tl, 'Real seeds: simulation vs tracked drop (validation only)');
    exportgraphics(f, fullfile(outDir, 'real_seeds_descent.png'), 'Resolution', 130);  close(f);

    % --- paths, side (height vs time) and top views, one tile per seed ---------
    sNums = unique([runs.seedNum]);
    for view = {'side', 'top'}
        f = figure('Visible', 'off', 'Position', [50 50 1500 1100]);
        tl = tiledlayout(f, 'flow', 'TileSpacing', 'tight');
        for s = sNums
            ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
            k = find(tNums == s, 1);
            if ~isempty(k)
                tr = tracks(k);
                if strcmp(view{1}, 'side'); plot(ax, tr.t, tr.pos(:,3), 'k.', 'MarkerSize', 3);
                else; plot(ax, tr.pos(:,1) - tr.pos(1,1), tr.pos(:,2) - tr.pos(1,2), 'k.', 'MarkerSize', 3); end
            end
            for i = 1:numel(rel)
                r = runs([runs.seedNum] == s & [runs.releaseDeg] == rel(i));
                if isempty(r) || isempty(r.pos); continue; end
                z0 = r.pos(1,3);  zt = z0;  if ~isempty(k); zt = tracks(k).pos(1,3); end
                if strcmp(view{1}, 'side'); plot(ax, r.tq, r.pos(:,3) - z0 + zt, '-', 'Color', col(i,:), 'LineWidth', 1);
                else; plot(ax, r.pos(:,1) - r.pos(1,1), r.pos(:,2) - r.pos(1,2), '-', 'Color', col(i,:), 'LineWidth', 1); end
            end
            if strcmp(view{1}, 'top'); axis(ax, 'equal'); else; ylim(ax, [0, 10.5]); end
            title(ax, sprintf('seed %d', s), 'FontSize', 8);  ax.FontSize = 7;
        end
        if strcmp(view{1}, 'side')
            title(tl, sprintf(['Height (m) vs time (s): tracked (black dots); simulated, released at the ' ...
                'track''s first point: %s'], relLegend(rel, col)));
        else
            title(tl, sprintf('Horizontal path from release (m): tracked (black); simulated: %s', relLegend(rel, col)));
        end
        exportgraphics(f, fullfile(outDir, sprintf('real_seeds_paths_%s.png', view{1})), 'Resolution', 110);
        close(f);
    end

    % --- simulated labels: seeds x releases -----------------------------------
    pal = seedModePalette();
    f = figure('Visible', 'off', 'Position', [100 100 520 900]);
    ax = axes(f);  hold(ax, 'on');
    for a = 1:numel(sNums)
        for i = 1:numel(rel)
            r = runs([runs.seedNum] == sNums(a) & [runs.releaseDeg] == rel(i));
            if isempty(r); continue; end
            j = find(strcmp(pal.names, char(r.simLabel)), 1);
            c = [0.6 0.6 0.6];  if ~isempty(j); c = pal.colors(j,:); end
            rectangle(ax, 'Position', [i-0.5, a-0.5, 1, 1], 'FaceColor', c, 'EdgeColor', 'w');
            text(ax, i, a, char(r.simLabel), 'HorizontalAlignment', 'center', 'FontSize', 7);
        end
    end
    set(ax, 'YDir', 'reverse', 'XTick', 1:numel(rel), 'XTickLabel', compose('%d deg', rel), ...
        'YTick', 1:numel(sNums), 'YTickLabel', compose('seed %d', sNums), 'FontSize', 8);
    xlim(ax, [0.5 numel(rel)+0.5]);  ylim(ax, [0.5 numel(sNums)+0.5]);
    title(ax, 'Simulated mode (six-part) by release');
    exportgraphics(f, fullfile(outDir, 'real_seeds_labels.png'), 'Resolution', 130);  close(f);
end

function s = relLegend(rel, col)
    names = {'blue', 'orange', 'yellow', 'purple', 'green', 'cyan', 'red'};
    s = strjoin(arrayfun(@(i) sprintf('%d deg = %s', rel(i), names{min(i, numel(names))}), ...
                1:numel(rel), 'UniformOutput', false), ', ');
end
