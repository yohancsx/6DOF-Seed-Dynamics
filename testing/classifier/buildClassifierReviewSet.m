%% Classifier review set -- runs, BLIND animations and a labelling spreadsheet
% Step 1 of the roadmap ("fix the classifier first"). This script's ONLY job is to
% produce the material for labelling flight modes by eye:
%
%   <outRoot>\<timestamp>_<githash>\
%       classifier_review.xlsx   sheet 'review'  : one row per run -- label_eye,
%                                                  confidence, notes, label_new
%                                                  (blank, for the new classifier),
%                                                  then the CURRENT code's labels
%                                sheet 'metrics' : every metric behind those labels
%                                sheet 'vocabulary': the eye labels allowed
%                                sheet 'about'   : config, git hash, how to use it
%       videos\R###.mp4          blind animation of the settled flight (see below)
%       summaries\R###.png       static summary: whole flight + attitude/rate traces
%       trajectories\R###.mat    the raw (t, x) and everything needed to rebuild the
%                                seed -- so a new classifier can be scored against
%                                these labels WITHOUT re-integrating, and the set
%                                can be frozen as a regression test later
%       manifest.mat
%
% WHAT IS RUN: a coarse chord x span nut-position grid on our seed (same seed,
% range and release as the baseline mode grid), under BOTH models, so a cell's
% planar and shape3d labels can be compared.
%
% BLIND BY DESIGN, so the eye labels are an independent check on the code:
%   * run IDs (R001 ...) are assigned in RANDOM order across both models and all
%     grid cells -- the ID gives away neither the model nor the nut position;
%   * videos use animateModeTrajectory's 'blind' option: no mode colours, no
%     legend, no label in any title;
%   * in the spreadsheet the eye columns come FIRST; the current code's labels sit
%     to the right of them. Hide columns H onward while labelling.
% Model and nut position are in the sheet (right-hand columns) for afterwards.
%
% Mode definitions, in words and math, and the eye vocabulary:
%   testing/classifier/MODE_DEFINITIONS.md   <- add new modes there.
%
% NEVER OVERWRITES: every run writes a new timestamped folder, so a spreadsheet
% you have already labelled cannot be clobbered by re-running this.
%
% COST: integration is cheap (parfor). The animations dominate -- about
% nRuns x animWindow x animFps frames. The defaults (72 runs x 2 s x 100 fps)
% are ~14,000 frames. Frame rate is in SIMULATED time, and it has to be high: an
% autorotating seed turns at ~10-16 rev/s, so at 20 fps its attitude aliases into
% nonsense. animPlaybackSpeed then slows the video down so it can be watched.
%
% Body axes: x = chord, y = normal, z = span. Section-organised, no local
% functions, so it saves as an .mlx cleanly.

%% 1. Configuration  -- EDIT HERE
root = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
outRoot = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Classifier Review";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% --- Which models, which grid ---------------------------------------------
cfg.models  = {'planar', 'shape3d'};   % each at its own builder defaults
cfg.gridN   = 6;                       % gridN x gridN nut positions per model
cfg.gridMax = 1.5;                     % nut offset range, x chord / x span (as baseline)

% --- Seed (our working seed, as in the baseline) --------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;
cfg.rhoFluid    = 1.225;   cfg.g = 9.81;

% --- Simulation + the classifier window -----------------------------------
cfg.tspan  = [0 12];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.q0     = axisAngleToQuat([0;0;1], pi/6);   % release, as the mode grids
cfg.omega0 = [0;0;0];
cfg.metricOpts = struct('windowStartFrac', 0.5, 'convergeTol', 0.20);
cfg.modeThresholds  = defaultModeThresholds();
cfg.paperThresholds = defaultPaperThresholds();

% --- Animations -----------------------------------------------------------
cfg.animate           = true;
cfg.animWindow        = 2;      % s of flight animated, counted back from the end
                                % (inside the classifier window, which is the last
                                %  (1-windowStartFrac)*tspan(2) = 6 s)
cfg.animFps           = 100;    % frames per SIMULATED second -- keep high, see header
cfg.animPlaybackSpeed = 0.15;   % video rate = animFps*this = 15 fps -> 13 s videos
cfg.randomSeed        = 1;      % fixes the run-ID shuffle, so a rerun is comparable

%% 2. Output folder, run list (shuffled, blind IDs)
[~, gh] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
gh = strtrim(gh);  if isempty(gh) || contains(gh, ' '); gh = 'nogit'; end
stamp  = char(datetime('now','Format','yyyy-MM-dd_HHmmss'));
outDir = fullfile(char(outRoot), sprintf('%s_%s', stamp, gh));
vidDir = fullfile(outDir, 'videos');  sumDir = fullfile(outDir, 'summaries');
trjDir = fullfile(outDir, 'trajectories');
cellfun(@mkdir, {outDir, vidDir, sumDir, trjDir});

xh = cfg.spanLength/2;   yh = cfg.chordLength/2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

frac = linspace(0, cfg.gridMax, cfg.gridN);
[CF, SF, MI] = ndgrid(frac, frac, 1:numel(cfg.models));   % chord, span, model
nRuns = numel(CF);
rng(cfg.randomSeed);
order = randperm(nRuns).';                % run k is the order(k)-th grid entry
CF = CF(:);  SF = SF(:);  MI = MI(:);     % flatten first: keeps every list nRuns x 1
chordFrac = CF(order);   spanFrac = SF(order);   modelIdx = MI(order);
runId = compose('R%03d', (1:nRuns).');
fprintf('\n=== classifier review set: %d runs -> %s\n', nRuns, outDir);

%% 3. Integrate every run (parfor)
tAll = cell(nRuns,1);  xAll = cell(nRuns,1);  status = repmat({'OK'}, nRuns, 1);
errMsg = repmat({''}, nRuns, 1);
odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);
models = cfg.models;
try if isempty(gcp('nocreate')); parpool('local'); end; catch; end
tInt = tic;
parfor k = 1:nRuns
    b = baseBsp;  b.tSamples = 0;  b.nutMass_t = cfg.nutMass;
    b.nutPos_t = [chordFrac(k)*cfg.chordLength; 0; spanFrac(k)*cfg.spanLength];
    c = cfg;  c.shapeModel = models{modelIdx(k)};
    try
        sp  = buildSeedParams(b, c);
        rhs = seedRHS(sp);
        x0  = [zeros(3,1); cfg.q0; zeros(3,1); cfg.omega0];
        [t, x] = ode45(@(tt,xx) rhs(tt,xx,sp), cfg.tspan, x0, odeOpts);
        tAll{k} = t;  xAll{k} = x;
        if any(~isfinite(x(:))); status{k} = 'NONFINITE'; end
    catch ME
        status{k} = 'FAILED';  errMsg{k} = ME.message;
    end
end
fprintf('  integrated %d runs in %.0f s (%d not OK)\n', nRuns, toc(tInt), ...
        sum(~strcmp(status,'OK')));

%% 4. Per run: classify (current code), save trajectory, summary figure, video
labCur  = strings(nRuns,1);  whyCur  = strings(nRuns,1);
labPap  = strings(nRuns,1);  whyPap  = strings(nRuns,1);
metricNames = {'descentSpeed','glideRatio','helixRadius','helixValid','coneAngleDeg', ...
               'tiltStd','verticalSpinMag','spanwiseSpin','tumbleFrac','converged', ...
               'spanAxisTiltDeg','spanAxisTiltStd','revolutionRate','selfRotationRate', ...
               'maxTurnsOneWay','tumbleReversals'};
M = nan(nRuns, numel(metricNames));
aopt = struct('window', cfg.animWindow, 'fps', cfg.animFps, ...
              'playbackSpeed', cfg.animPlaybackSpeed, 'showSeedVels', false, ...
              'blind', true, 'topPanel', 'groundTrack', 'savePng', false);
tAnim = tic;
for k = 1:nRuns
    mdl = models{modelIdx(k)};
    nut = [chordFrac(k)*cfg.chordLength; 0; spanFrac(k)*cfg.spanLength];
    t = tAll{k};  x = xAll{k};

    % --- current classifiers (a nonfinite/failed run gets the grid's label) --
    if strcmp(status{k}, 'OK')
        pm = computePaperMetrics(t, x, cfg.metricOpts);
        [lc, ic] = classifyFlightMode(pm, cfg.modeThresholds);
        [lp, ip] = classifyPaperMode(pm, cfg.paperThresholds);
        labCur(k) = lc;  whyCur(k) = ic.reason;  labPap(k) = lp;  whyPap(k) = ip.reason;
        for j = 1:numel(metricNames);  M(k,j) = double(pm.(metricNames{j}));  end
    elseif strcmp(status{k}, 'NONFINITE')
        labCur(k) = "chaotic";  whyCur(k) = "non-finite state (grid convention)";
        labPap(k) = "failed";   whyPap(k) = "non-finite state";
    else
        labCur(k) = "failed";   whyCur(k) = string(errMsg{k});
        labPap(k) = "failed";   whyPap(k) = string(errMsg{k});
    end

    % --- raw trajectory + what is needed to rebuild the seed --------------
    runInfo = struct('runId', runId{k}, 'model', mdl, 'chordFrac', chordFrac(k), ...
                     'spanFrac', spanFrac(k), 'nutPos', nut, 'status', status{k}, ...
                     'labelCurrent', labCur(k), 'labelPaperCurrent', labPap(k)); %#ok<NASGU>
    save(fullfile(trjDir, [runId{k} '.mat']), 't', 'x', 'runInfo', 'cfg', 'baseBsp');

    if ~strcmp(status{k}, 'OK')
        fprintf('  [%3d/%d] %s  %-7s  %s\n', k, nRuns, runId{k}, mdl, status{k});
        continue
    end

    % --- static summary: whole flight + attitude and rate traces ----------
    % Shows what a short video cannot: the whole descent, and whether the
    % attitude is steady, rocking or turning over across the full window.
    nS = numel(t);  th = zeros(nS,1);  ga = zeros(nS,1);  wV = zeros(nS,1);
    for i = 1:nS
        R = quatToRotm(x(i,4:7).' / norm(x(i,4:7)));
        s = R*[0;0;1];  n = R*[0;1;0];  wW = R*x(i,11:13).';
        th(i) = asind(max(-1,min(1,s(2))));   ga(i) = acosd(min(1,abs(n(2))));
        wV(i) = wW(2);
    end
    inW = t >= cfg.metricOpts.windowStartFrac * t(end);
    f = figure('Visible','off','Color','w','Position',[50 50 1800 800]);
    tl = tiledlayout(f, 2, 4, 'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('%s   (blind: model and nut position are in the spreadsheet)', ...
                      runId{k}), 'FontWeight','bold');
    ax = nexttile(tl, 1, [2 1]);
    plot3(ax, x(:,1), x(:,3), x(:,2), '-', 'Color', [0.75 0.75 0.75]);  hold(ax,'on');
    plot3(ax, x(inW,1), x(inW,3), x(inW,2), '-', 'Color', [0.15 0.35 0.75], 'LineWidth', 1.5);
    grid(ax,'on');  axis(ax,'equal');  view(ax, [40 18]);
    xlabel(ax,'X (m)'); ylabel(ax,'Z (m)'); zlabel(ax,'Y up (m)');
    title(ax, 'whole flight (blue = classifier window)');
    ax = nexttile(tl, 2);
    plot(ax, x(inW,1), x(inW,3), 'Color', [0.15 0.35 0.75]);  axis(ax,'equal');  grid(ax,'on');
    xlabel(ax,'X (m)'); ylabel(ax,'Z (m)');  title(ax, 'ground track, classifier window');
    ax = nexttile(tl, 3);
    plot(ax, t, th, t, ga);  grid(ax,'on');  xline(ax, t(find(inW,1)), 'k:');
    ylim(ax, [-90 90]);  legend(ax, {'\theta span axis from horizontal', ...
         '\gamma normal from vertical'}, 'Location','best');
    xlabel(ax,'t (s)');  ylabel(ax,'deg');  title(ax, 'attitude');
    % Cumulative rotation, in turns. The clearest single tell between the
    % rotating modes: tumbling climbs steadily about the span axis, fluttering
    % stays bounded, revolution climbs about the vertical. (Raw integral of
    % omega_z: for a revolving seed it includes the revolution's projection on
    % the tilted span axis -- see MODE_DEFINITIONS section 1.)
    ax = nexttile(tl, 4);
    plot(ax, t, cumtrapz(t, x(:,13))/(2*pi), t, cumtrapz(t, wV)/(2*pi));
    grid(ax,'on');  xline(ax, t(find(inW,1)), 'k:');
    legend(ax, {'about span axis (\int\omega_z)','about world vertical'}, 'Location','best');
    xlabel(ax,'t (s)');  ylabel(ax,'turns');  title(ax, 'cumulative rotation');
    ax = nexttile(tl, 6);
    plot(ax, t, x(:,11:13));  grid(ax,'on');  xline(ax, t(find(inW,1)), 'k:');
    legend(ax, {'\omega_x chord','\omega_y normal','\omega_z span'}, 'Location','best');
    xlabel(ax,'t (s)');  ylabel(ax,'rad/s');  title(ax, 'body rates');
    ax = nexttile(tl, 7);
    plot(ax, t, wV, t, -x(:,9));  grid(ax,'on');  xline(ax, t(find(inW,1)), 'k:');
    legend(ax, {'spin about world vertical (rad/s)','descent speed (m/s)'}, 'Location','best');
    xlabel(ax,'t (s)');  title(ax, 'revolution and descent');
    exportgraphics(f, fullfile(sumDir, [runId{k} '.png']), 'Resolution', 110);
    close(f);

    % --- blind animation of the settled flight ----------------------------
    if cfg.animate
        c = cfg;  c.shapeModel = mdl;
        b = baseBsp;  b.tSamples = 0;  b.nutMass_t = cfg.nutMass;  b.nutPos_t = nut;
        sp = buildSeedParams(b, c);
        animateCaseVideo(t, x, sp, runId{k}, vidDir, aopt);
    end
    fprintf('  [%3d/%d] %s done  (%.0f s elapsed)\n', k, nRuns, runId{k}, toc(tAnim));
end

%% 5. The spreadsheet
video   = strcat("videos\", string(runId), ".mp4");
summary = strcat("summaries\", string(runId), ".png");
blank   = strings(nRuns,1);
review = table(string(runId), video, summary, blank, nan(nRuns,1), blank, blank, ...
               labCur, labPap, whyCur, whyPap, string(models(modelIdx)).', ...
               chordFrac, spanFrac, 1e3*chordFrac*cfg.chordLength, ...
               1e3*spanFrac*cfg.spanLength, string(status), ...
    'VariableNames', {'run_id','video','summary','label_eye','confidence_eye_1to3', ...
                      'notes_eye','label_new','label_current','label_paper_current', ...
                      'reason_current','reason_paper','model','nut_chordFrac', ...
                      'nut_spanFrac','nut_x_mm','nut_z_mm','status'});
M(~isfinite(M)) = NaN;       % Inf (e.g. no helix) would land in Excel as 65535
metrics = [table(string(runId), 'VariableNames', {'run_id'}), ...
           array2table(M, 'VariableNames', metricNames)];
vocab = table( ...
    ["gliding";"diving";"parachuting";"falling";"fluttering";"tumbling"; ...
     "spiralTumbling";"autorotation";"chaotic";"endOn";"edgeSlide";"unsure"], ...
    ["no rotation, sideways travel >= height lost"; ...
     "no rotation, fast steep fall, chord edge down"; ...
     "no rotation, broadside, slow"; ...
     "no rotation, steady large tilt heavy side down (paper FA)"; ...
     "rocks about the span axis, never a full turn (paper FL)"; ...
     "whole end-over-end turns about the span, no revolution"; ...
     "tumbles while spiralling about a vertical axis (paper CST/SST)"; ...
     "revolves about a vertical axis at a steady attitude (paper AR)"; ...
     "no sustained pattern (paper CH)"; ...
     "NON-PHYSICAL: span axis near vertical, falls span-first"; ...
     "NON-PHYSICAL: slides along its own span, edge-on"; ...
     "cannot tell -- say why in notes_eye"], ...
    'VariableNames', {'label','meaning'});
about = table(["created"; "git"; "models"; "grid"; "release"; "tspan"; ...
               "classifier window"; "videos"; "definitions"; "how to use"], ...
    [string(stamp); string(gh); strjoin(string(cfg.models), ', '); ...
     sprintf('%dx%d, nut 0..%.1f chord and span', cfg.gridN, cfg.gridN, cfg.gridMax); ...
     "pi/6 about body z, from rest"; sprintf('%g s', cfg.tspan(2)); ...
     sprintf('t >= %.2f t_end', cfg.metricOpts.windowStartFrac); ...
     sprintf('last %g s, %g fps sim time, played at %gx', cfg.animWindow, cfg.animFps, ...
             cfg.animPlaybackSpeed); ...
     "testing/classifier/MODE_DEFINITIONS.md (eye vocabulary: section 4)"; ...
     "Hide columns H onward (current labels), watch each video + summary, fill label_eye / confidence / notes. label_new is filled later by the new classifier."], ...
    'VariableNames', {'item','value'});
xls = fullfile(outDir, 'classifier_review.xlsx');
writetable(review,  xls, 'Sheet', 'review');
writetable(metrics, xls, 'Sheet', 'metrics');
writetable(vocab,   xls, 'Sheet', 'vocabulary');
writetable(about,   xls, 'Sheet', 'about');
save(fullfile(outDir, 'manifest.mat'), 'cfg', 'stamp', 'gh', 'runId', 'chordFrac', ...
     'spanFrac', 'modelIdx', 'status', 'labCur', 'labPap', 'metricNames', 'M');

% --- quick look at what the current code says ---------------------------
fprintf('\nCurrent-code labels (for reference -- do not look before labelling):\n');
for im = 1:numel(models)
    sel = (modelIdx == im) & strcmp(status, 'OK');     % both nRuns x 1
    if ~any(sel); fprintf('  %-8s (no finite runs)\n', models{im}); continue; end
    [u,~,j] = unique(labCur(sel));
    fprintf('  %-8s %s\n', models{im}, strjoin(compose('%s=%d', u, accumarray(j,1)), '  '));
end
fprintf('\nWritten to %s\n', outDir);
