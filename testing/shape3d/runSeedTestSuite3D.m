%% Seed dynamics test suite (3D MODEL) -- one-click sweep runner
% The shape3d counterpart of testing/planar/runSeedTestSuite.m, and called the
% same way: edit the configuration block below and run the file. Everything it
% writes carries a 3D marker, so 2D and 3D outputs can never be confused.
%
% WHAT IT RUNS (each stage can be switched off in section 1):
%   A. TWIST sweep      -- spanwise geometric pitch, 0..35 deg tip, centred nut,
%                          released LEVEL (cfg.q0Twist) so the spin-up is the
%                          shape's doing and not the release tilt's.
%                          One ANIMATION per case.
%   B. CURVATURE sweep  -- spanwise dihedral, two families, centred nut, released
%                          at cfg.q0. One ANIMATION per case.
%   C. SWEEP SUITE      -- the same nine groups the 2D suite runs (nut mass, nut
%                          chord/span/diagonal position, initial pitch, initial
%                          roll, initial yaw spin, one-sided lift/drag asymmetry,
%                          strip-count convergence), through the SAME
%                          seedTestCases / runOneSeedCase helpers, only with
%                          shapeModel = 'shape3d'. Per-case .fig/.png/.mat and a
%                          per-group overlay, exactly as in 2D, so a group can be
%                          laid next to its planar twin.
%   D. CoM MOVEMENT     -- the three scripted moving-nut scenarios, exercising the
%                          time-varying mass path (nutPos_t -> com_t, I_G_t,
%                          I_G_dot_t) under the 3D model. One ANIMATION each.
%   E. MODE GRID        -- chord x span nut-position map on OUR seed, 3D model.
%                          NOTE this is the grid that does not otherwise exist as
%                          a standalone script: testing/planar/runSeedModeGrid.m
%                          calls seed6DOFODE directly and is planar-only, so it
%                          contains NONE of the 3D terms.
%   F. PAPER COMPARISON -- ON BY DEFAULT. Hou et al. (2025): their plate, their
%                          parameter window, their release, the five published
%                          anchor cases and the phase map, scored against the
%                          digitised Fig. 2a.
%
% MODE LABELS: the six-part classifier (classifySeedMode; definitions in
% testing/classifier/MODE_DEFINITIONS.md) is the label of record in every stage;
% the older classifyFlightMode label is kept beside it in [brackets].
%
% TRACKED NUMBERS (stages E and F, via summarizeModeGrid): paper agreement,
% non-physical fraction, and the span-tilt distribution of revolvers. Written to
% <run>/tracking_3D.txt, and one row per run is appended to
% <outPath>/tracking_history_3D.csv so runs can be compared over time.
% Run folders are named <timestamp>_<githash>_3D_<label>.
%
% WHY IT EXISTS: to turn physics terms off and see what happens. Every stage runs
% with the SAME switch overrides, applied through cfg so buildSeedParams warns if
% a switch does not exist in this model rather than silently ignoring it. Give
% each configuration a label and the run folders can be diffed directly.
%
% COST: the animations dominate the wall clock, not the integrations. Stages A/B/D
% render ~35 videos; cfg.animWindow keeps each one to the last few seconds of the
% flight (where the settled behaviour is) rather than the whole drop. Set
% cfg.animate = false for a numbers-only run.
%
% Body axes: x = chord, y = normal, z = span.
% Section-organised and free of local functions, so it saves as an .mlx cleanly.

%% 1. Configuration  -- EDIT HERE
root    = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
outPath = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\test_suite_3D";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% --- Label this configuration (goes in the run folder name) ---------------
cfg.label = 'all_enabled';           % e.g. 'noKirchhoff', 'LEVon', 'noEdgeDrag'

% --- PHYSICS SWITCHES -- the point of this suite --------------------------
% Set a field to override the builder default; COMMENT IT OUT to use the default.
% Whatever is left at default is recorded in the manifest, read off a built seed,
% so a run always says which physics it actually used.
sw = struct();
% sw.enableAddedMassRate   = false;   % Adot*v      (APW eqs. 6.1-6.2 cross terms)
% sw.enableAddedMassMoment = false;   % v x (A v)   (APW eq. 6.3, the Munk moment)
% sw.enableEdgeDrag        = false;   % tip crossflow drag
sw.enableLEV             = true;    % leading-edge-vortex vortex lift
sw.levApplicationPoint   = 'forward';
% sw.lev = struct('rossbyDefinition','geometric');

% --- Which stages to run --------------------------------------------------
cfg.runTwist      = true;
cfg.runCurvature  = true;
cfg.runSweeps     = true;      % the nine 2D-parity groups
cfg.runComMove    = true;      % moving-nut scenarios
cfg.runModeGrid   = true;
cfg.runPaperGrid  = true;      % ON by default -- the comparison that matters

% --- Animations (stages A, B, D) ------------------------------------------
cfg.animate           = true;
cfg.animWindow        = Inf;     % seconds animated, counted back from the END (Inf = all)
cfg.animFps           = 20;
cfg.animPlaybackSpeed = 0.25;  % <1 = slow motion, so a fast spin is actually visible
cfg.animShowVels      = false; % per-strip wind arrows in the zoom panel (slow)

% --- Baseline seed geometry / material (our working seed) -----------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   cfg.nutPos      = [0;0;0];
cfg.rhoFluid    = 1.225;   cfg.g = 9.81;
cfg.shapeModel  = 'shape3d';

% --- Simulation -----------------------------------------------------------
cfg.tspan = [0 12];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();          % older classifier (continuity)
cfg.seedModeThresholds = defaultSeedModeThresholds();  % six-part classifier (of record)
% Every mode-grid run appends one row of tracked numbers (paper agreement,
% non-physical fraction, revolver span tilt) to this file, so runs can be compared
% over time:
cfg.trackingFile = fullfile(char(outPath), 'tracking_history_3D.csv');

% --- Shape sweep settings (stages A, B) -----------------------------------
cfg.tipDeg    = [0 5 10 15 20 25 30 35];      % twist / curvature tip angle (deg)
cfg.q0        = axisAngleToQuat([0;0;1], pi/6);  % release for curvature, CoM move, our grid
% TWIST releases LEVEL. A twisted seed already carries a built-in geometric angle
% of attack that varies along the span, so adding a pitched release confounds the
% two: the spin-up would be partly the shape and partly the initial tilt, and the
% 0 deg row would not be the flat reference it is meant to be.
cfg.q0Twist   = [1;0;0;0];
cfg.omega0    = [0;0;0];

% --- Sweep-suite settings (stage C) -- same values as the 2D suite ---------
cfg.sweepTspan    = [0 5];      % 2D suite uses 5 s; kept so groups compare directly
cfg.nIncr         = 5;          % increments per numeric sweep
cfg.nutMassFrac   = 0.20;       % nut mass +/- 20%
cfg.nutPosMaxFrac = 1.20;       % nut position 0 .. 120% of half-dimension
cfg.tiltMaxDeg    = 45;         % pitch / roll release tilt, -45 .. +45 deg
cfg.yawSpinMin    = 1;          % initial spin about the plate normal (rad/s)
cfg.yawSpinMax    = 5;
cfg.asymFactor    = 0.5;        % one-sided lift/drag multiplier
cfg.stripCounts   = [1 5 10 20];% strip-convergence counts
cfg.groups = struct('nutMass',true, 'nutChord',true, 'nutSpan',true, 'nutDiag',true, ...
                    'pitch',true, 'roll',true, 'yawSpin',true, 'asymmetry',true, ...
                    'stripConv',true);
cfg.savePng = true;   cfg.saveOverlay = true;

% --- CoM movement settings (stage D) --------------------------------------
cfg.comDwell    = 1.0;    % pause held at each nut waypoint (s)
cfg.comMove     = 0.8;    % transition time between waypoints (s)
cfg.comDt       = 0.02;   % dense nut-path / mass-sample step (s)
cfg.comSpanMax  = 1.0;    % scenario 1 amplitude, x half-span
cfg.comChordMax = 1.0;    % scenario 2 amplitude, x half-chord

% --- Mode grid + paper grid (stages E, F) ---------------------------------
cfg.gridN     = 20;                            % mode grid is gridN x gridN
cfg.gridMax   = 1.5;                           % nut offset range, in chord/span units
cfg.paperNX   = 13;   cfg.paperNY = 8;         % paper-window grid resolution

%% 2. Apply the switch overrides, build the run folder
fn = fieldnames(sw);
for k = 1:numel(fn);  cfg.(fn{k}) = sw.(fn{k});  end

timestamp = char(datetime('now','Format','yyyy-MM-dd_HHmmss'));
[gst, gh] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
gh = strtrim(gh);  if gst ~= 0 || isempty(gh) || contains(gh, ' '); gh = 'nogit'; end
[~, gd] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', root));
if ~isempty(strtrim(gd)); gh = [gh '-dirty']; end
% <timestamp>_<githash>_3D_<label>: the hash pins the code, as for baseline snapshots
runDir = fullfile(char(outPath), sprintf('%s_%s_3D_%s', timestamp, gh, cfg.label));
mkdir(runDir);

% Mode labels: the six-part classifier (classifySeedMode) is the label of record;
% the older classifyFlightMode label is kept beside it for continuity with past runs.
so = cfg.metricOpts;   so.spanLength = cfg.spanLength;
so.flipHysteresisTurns = cfg.seedModeThresholds.flipHysteresisTurns;

xh = cfg.spanLength/2;   yh = cfg.chordLength/2;   hs = xh;   hc = yh;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

aopt = struct('window',cfg.animWindow, 'fps',cfg.animFps, ...
              'playbackSpeed',cfg.animPlaybackSpeed, 'showSeedVels',cfg.animShowVels);

% Record the physics ACTUALLY in force, read off a built seed rather than assumed.
bspProbe = baseBsp;  bspProbe.tSamples = 0;
bspProbe.nutPos_t = [0;0;0];  bspProbe.nutMass_t = cfg.nutMass;
spProbe = buildSeedParams(bspProbe, cfg);
effSwitches = struct('enableEdgeDrag', spProbe.enableEdgeDrag, ...
                     'enableAddedMassRate', spProbe.enableAddedMassRate, ...
                     'enableAddedMassMoment', spProbe.enableAddedMassMoment, ...
                     'enableLEV', spProbe.enableLEV, ...
                     'levApplicationPoint', spProbe.levApplicationPoint, ...
                     'rossbyDefinition', spProbe.lev.rossbyDefinition);
fprintf('\n=== 3D test suite: %s ===\n', runDir);
fprintf('physics in force: edgeDrag %d, addedMassRate %d, addedMassMoment %d, LEV %d (%s, %s)\n', ...
        effSwitches.enableEdgeDrag, effSwitches.enableAddedMassRate, ...
        effSwitches.enableAddedMassMoment, effSwitches.enableLEV, ...
        effSwitches.levApplicationPoint, effSwitches.rossbyDefinition);

summary = strings(0,1);
summary(end+1) = sprintf('3D test suite  %s  git %s  label=%s', timestamp, gh, cfg.label); %#ok<*SAGROW>
summary(end+1) = "mode labels: six-part classifier (classifySeedMode); older classifyFlightMode label in [brackets]";
summary(end+1) = sprintf('physics: edgeDrag %d, rate %d, moment %d, LEV %d (%s, %s)', ...
        effSwitches.enableEdgeDrag, effSwitches.enableAddedMassRate, ...
        effSwitches.enableAddedMassMoment, effSwitches.enableLEV, ...
        effSwitches.levApplicationPoint, effSwitches.rossbyDefinition);

%% 3. Stage A -- TWIST sweep (+ one animation per case)
twist = struct('tipDeg',[], 'mode',strings(0,1), 'descent',[], 'vSpin',[], 'cone',[]);
if cfg.runTwist
    fprintf('\n--- A. twist sweep ---\n');
    twDir = fullfile(runDir, 'twist_3D');   mkdir(twDir);
    nT = numel(cfg.tipDeg);
    tw_mode = strings(nT,1);  tw_desc = nan(nT,1);  tw_vs = nan(nT,1);  tw_cone = nan(nT,1);
    tw_new = strings(nT,1);
    for j = 1:nT
        kt = deg2rad(cfg.tipDeg(j))/hs;
        bb = baseBsp;  bb.twist = @(z) kt .* z;      % anti-symmetric twist
        nm = sprintf('twist_%02ddeg', cfg.tipDeg(j));
        r = runSingleMode(nm, [0;0;0], cfg.q0Twist, cfg.omega0, cfg, bb);
        tw_new(j)  = string(classifySeedMode(computeSeedModeMetrics(r.t, r.x, so), ...
                                             cfg.seedModeThresholds));
        tw_mode(j) = string(r.mode);   tw_desc(j) = r.metrics.descentSpeed;
        tw_vs(j)   = r.metrics.verticalSpinMag;  tw_cone(j) = r.metrics.coneAngleDeg;
        if cfg.animate
            animateCaseVideo(r.t, r.x, r.seedParams, nm, twDir, aopt);
        end
    end
    twist = struct('tipDeg',cfg.tipDeg, 'mode',tw_new, 'modeOld',tw_mode, ...
                   'descent',tw_desc, 'vSpin',tw_vs, 'cone',tw_cone);
    save(fullfile(runDir,'twist_sweep_3D.mat'), 'twist', 'cfg', 'effSwitches');
    summary(end+1) = "";
    summary(end+1) = "TWIST sweep (tip deg -> mode [old], descent, vSpin, cone)";
    for j = 1:nT
        summary(end+1) = sprintf('  %2d  %-18s [%-12s] %6.2f %7.2f %6.1f', cfg.tipDeg(j), ...
                                 tw_new(j), tw_mode(j), tw_desc(j), tw_vs(j), tw_cone(j));
    end
end

%% 4. Stage B -- CURVATURE sweep, two families (+ one animation per case)
curv = struct();
if cfg.runCurvature
    fprintf('\n--- B. curvature sweep ---\n');
    cvDir = fullfile(runDir, 'curvature_3D');   mkdir(cvDir);
    nC = numel(cfg.tipDeg);
    cv_mode = strings(nC,2);  cv_desc = nan(nC,2);  cv_vs = nan(nC,2);  cv_comY = nan(nC,2);
    cv_new  = strings(nC,2);
    famName = {'bowl','oneSide'};
    for fam = 1:2
        for j = 1:nC
            kc = deg2rad(cfg.tipDeg(j))/hs;
            bb = baseBsp;
            if fam == 1
                bb.curvature = @(s) kc .* s;                 % symmetric bowl
            else
                bb.curvature = @(s) (s > 0) .* kc .* s;      % one side only
            end
            nm = sprintf('%s_%02ddeg', famName{fam}, cfg.tipDeg(j));
            r = runSingleMode(nm, [0;0;0], cfg.q0, cfg.omega0, cfg, bb);
            cv_new(j,fam)  = string(classifySeedMode(computeSeedModeMetrics(r.t, r.x, so), ...
                                                     cfg.seedModeThresholds));
            cv_mode(j,fam) = string(r.mode);   cv_desc(j,fam) = r.metrics.descentSpeed;
            cv_vs(j,fam)   = r.metrics.verticalSpinMag;
            cv_comY(j,fam) = r.seedParams.massParams.com_t(2,1);
            if cfg.animate
                animateCaseVideo(r.t, r.x, r.seedParams, nm, cvDir, aopt);
            end
        end
    end
    curv = struct('tipDeg',cfg.tipDeg, 'family',{famName}, 'mode',cv_new, ...
                  'modeOld',cv_mode, 'descent',cv_desc, 'vSpin',cv_vs, 'comY',cv_comY);
    save(fullfile(runDir,'curvature_sweep_3D.mat'), 'curv', 'cfg', 'effSwitches');
    summary(end+1) = "";
    summary(end+1) = "CURVATURE sweep (tip deg -> bowl mode [old] descent | oneSide mode [old] descent)";
    for j = 1:nC
        summary(end+1) = sprintf('  %2d  %-18s [%-12s] %6.2f | %-18s [%-12s] %6.2f', ...
            cfg.tipDeg(j), cv_new(j,1), cv_mode(j,1), cv_desc(j,1), ...
            cv_new(j,2), cv_mode(j,2), cv_desc(j,2));
    end
end

%% 5. Stage C -- SWEEP SUITE: the nine groups the 2D suite runs
% Built from the SAME seedTestCases / runOneSeedCase helpers as
% testing/planar/runSeedTestSuite.m -- the only difference is cfg.shapeModel, so
% any divergence from the planar run is the 3D physics and nothing else. Output
% layout mirrors the 2D suite (per-group subfolder, per-case .mat/.fig/.png, a
% group overlay) so the two runs can be opened side by side.
if cfg.runSweeps
    fprintf('\n--- C. sweep suite (2D-parity groups) ---\n');
    swDir = fullfile(runDir, 'sweeps_3D');   mkdir(swDir);
    scfg = cfg;   scfg.tspan = cfg.sweepTspan;
    [cases, swBsp] = seedTestCases(scfg);
    swBaseSp = buildSeedParams(swBsp, scfg);
    fprintf('  %d cases across enabled groups\n', numel(cases));
    swResults = [];
    for k = 1:numel(cases)
        cs = cases(k);
        fprintf('  [%2d/%2d] %-10s %-24s ... ', k, numel(cases), cs.group, cs.label);
        r = runOneSeedCase(cs, swBaseSp, swBsp, scfg);
        fprintf('%s\n', r.status);
        caseDir = fullfile(swDir, r.group);
        if ~exist(caseDir,'dir'); mkdir(caseDir); end
        result = r; %#ok<NASGU>
        save(fullfile(caseDir, [r.label '.mat']), 'result');
        if ~isempty(r.x)
            fig = plotSeedCaseTrajectory(r, 'off');
            savefig(fig, fullfile(caseDir, [r.label '.fig']));
            if cfg.savePng
                exportgraphics(fig, fullfile(caseDir, [r.label '.png']), 'Resolution',150);
            end
            close(fig);
        end
        swResults = [swResults; r]; %#ok<AGROW>
    end
    if cfg.saveOverlay
        gNames = unique({swResults.group}, 'stable');
        for gi = 1:numel(gNames)
            gRes = swResults(strcmp({swResults.group}, gNames{gi}));
            fig  = plotGroupOverlay(gRes, gNames{gi}, 'off');
            savefig(fig, fullfile(swDir, gNames{gi}, '_overlay.fig'));
            if cfg.savePng
                exportgraphics(fig, fullfile(swDir, gNames{gi}, '_overlay.png'), 'Resolution',150);
            end
            close(fig);
        end
    end
    sweepSummary = struct('group',{swResults.group}, 'label',{swResults.label}, ...
                          'status',{swResults.status}, 'sweepVal',{swResults.sweepVal}, ...
                          'descended',{swResults.descended});
    save(fullfile(runDir,'sweep_suite_3D.mat'), 'sweepSummary', 'scfg', 'effSwitches');
    nBad = sum(~strcmp({swResults.status}, 'OK'));
    fprintf('  %d cases, %d not OK\n', numel(swResults), nBad);
    summary(end+1) = "";
    summary(end+1) = sprintf('SWEEP SUITE -- %d cases, %d not OK', numel(swResults), nBad);
    summary(end+1) = sprintf('  %-11s %-24s %-10s %s', 'group','label','status','descended');
    for k = 1:numel(swResults)
        summary(end+1) = sprintf('  %-11s %-24s %-10s %d', swResults(k).group, ...
                                 swResults(k).label, swResults(k).status, swResults(k).descended);
    end
end

%% 6. Stage D -- CoM MOVEMENT scenarios (+ one animation each)
% The only stage that exercises the TIME-VARYING mass path under the 3D model:
% a dense nutPos_t schedule drives com_t / I_G_t / I_G_dot_t, and the I_G_dot*w
% term in rotationDynamics goes live. Each dwell is classified once the seed has
% settled into it (the last 55% of the dwell), so the scenario reads as a
% sequence of modes rather than one averaged blur.
%
% CAVEAT (unchanged from the planar version): sliding an internal mass on a
% free-flying body is not momentum-rigorous here -- the model prescribes where the
% CoM sits in the body, but adds no reaction from the internal mass's motion.
if cfg.runComMove
    fprintf('\n--- D. CoM movement ---\n');
    cmDir = fullfile(runDir, 'com_movement_3D');   mkdir(cmDir);
    clear scen
    scen(1) = struct('name','S1_spanwise_sweep', ...
        'pos',{{[0;0;0],[0;0;+cfg.comSpanMax*hs],[0;0;0],[0;0;-cfg.comSpanMax*hs],[0;0;0]}});
    scen(2) = struct('name','S2_chordwise_sweep', ...
        'pos',{{[0;0;0],[+cfg.comChordMax*hc;0;0],[0;0;0],[-cfg.comChordMax*hc;0;0],[0;0;0]}});
    scen(3) = struct('name','S3_glide_spin_dive', ...
        'pos',{{[0.8*hc;0;0],[0;0;0],[0;0;1.2*hs],[1.5*hc;0;0]}});
    cmOde = odeset('RelTol',cfg.odeRelTol,'AbsTol',cfg.odeAbsTol);
    summary(end+1) = "";
    summary(end+1) = "CoM MOVEMENT (scenario -> settled mode at each nut waypoint)";
    for si = 1:numel(scen)
        posList = [scen(si).pos{:}];
        dwell   = cfg.comDwell * ones(1, size(posList,2));
        [tD, nutPos, eventT, dwellInt] = buildComPath(posList, dwell, cfg.comMove, cfg.comDt);
        bb = baseBsp;  bb.tSamples = tD;  bb.nutPos_t = nutPos;
        bb.nutMass_t = cfg.nutMass * ones(size(tD));
        sp  = buildSeedParams(bb, cfg);
        rhs = seedRHS(sp);
        x0  = [zeros(3,1); cfg.q0; zeros(3,1); cfg.omega0];
        [t, x] = ode45(@(tt,xx) rhs(tt,xx,sp), [tD(1) tD(end)], x0, cmOde);
        dwellModes = strings(size(dwellInt,1),1);  dwellModesOld = dwellModes;
        soD = so;  soD.windowStartFrac = 0;          % the dwell slice IS the window
        for i = 1:size(dwellInt,1)
            sel = t >= dwellInt(i,1)+0.45*(dwellInt(i,2)-dwellInt(i,1)) & t <= dwellInt(i,2);
            if nnz(sel) >= 8
                dwellModes(i) = string(classifySeedMode( ...
                    computeSeedModeMetrics(t(sel), x(sel,:), soD), cfg.seedModeThresholds));
                m = computeTrajectoryMetrics(t(sel), x(sel,:), struct('windowStartFrac',0));
                % A dwell is short by construction, so do not also penalise it
                % for failing the convergence test that length would fail.
                m.converged = true;
                dwellModesOld(i) = string(classifyFlightMode(m, cfg.modeThresholds));
            else
                dwellModes(i) = "(short)";  dwellModesOld(i) = "(short)";
            end
        end
        save(fullfile(cmDir,[scen(si).name '.mat']), 't','x','tD','nutPos','eventT', ...
             'dwellInt','dwellModes','dwellModesOld','cfg','effSwitches');
        fprintf('  %-20s %.1f s, final height %+.3f m, modes: %s\n', scen(si).name, ...
                tD(end), x(end,2), strjoin(cellstr(dwellModes),' -> '));
        summary(end+1) = sprintf('  %-20s %s', scen(si).name, strjoin(cellstr(dwellModes),' -> '));
        summary(end+1) = sprintf('  %-20s [%s]', '', strjoin(cellstr(dwellModesOld),' -> '));
        if cfg.animate
            ao = aopt;  ao.window = Inf;  ao.eventTimes = eventT;   % the whole scripted motion
            animateCaseVideo(t, x, sp, scen(si).name, cmDir, ao);
        end
    end
end

%% 7. Stage E -- MODE GRID on our seed (3D model), six-part classifier
% Labels of record: classifySeedMode (testing/classifier/MODE_DEFINITIONS.md). The
% older classifyFlightMode label is saved beside it (modeIdxOld) so maps from
% earlier runs stay comparable. Tracked numbers (summarizeModeGrid): non-physical
% fraction and the span tilt of revolvers.
gridTrack = [];
if cfg.runModeGrid
    fprintf('\n--- E. mode grid (our seed, 3D) ---\n');
    chordFrac = linspace(0, cfg.gridMax, cfg.gridN);
    spanFrac  = linspace(0, cfg.gridMax, cfg.gridN);
    c = cfg.chordLength;  S = cfg.spanLength;
    pal = seedModePalette();
    modeList = pal.names;   modeCol = pal.colors;
    [modeListOld, ~] = seedModeColors();
    th = cfg.modeThresholds;   ths = cfg.seedModeThresholds;
    odeOpts = odeset('RelTol',cfg.odeRelTol,'AbsTol',cfg.odeAbsTol);
    total = cfg.gridN^2;  nXg = cfg.gridN;
    gLab = repmat({'failed'}, total, 1);  gOld = gLab;  gParts = repmat({''}, total, 1);
    gTilt = nan(total,1);  gRev = false(total,1);
    gDesc = nan(total,1);  gVs = nan(total,1);  gCone = nan(total,1);
    try if isempty(gcp('nocreate')); parpool('local'); end; catch; end
    tGrid = tic;
    parfor k = 1:total
        iz = floor((k-1)/nXg)+1;   ix = mod(k-1,nXg)+1;
        b = baseBsp;  b.tSamples = 0;
        b.nutPos_t = [chordFrac(ix)*c; 0; spanFrac(iz)*S];  b.nutMass_t = cfg.nutMass;
        lbl = 'failed';  old = 'failed';  prt = '';  tilt = NaN;  rev = false;
        d = NaN;  v = NaN;  cn = NaN;
        try
            sp  = buildSeedParams(b, cfg);
            rhs = seedRHS(sp);
            x0  = [zeros(3,1); cfg.q0; zeros(3,1); cfg.omega0];
            [t, x] = ode45(@(tt,xx) rhs(tt,xx,sp), cfg.tspan, x0, odeOpts);
            if all(isfinite(x(:)))
                sm = computeSeedModeMetrics(t, x, so);
                [lbl, p] = classifySeedMode(sm, ths);
                prt  = sprintf('%s/%s/%s/%s/%s', p.flip, p.revolution, p.attitude, ...
                               p.path, p.regularity);
                tilt = sm.absSpanTiltDeg;   rev = ~strcmp(p.revolution, 'none');
                old  = classifyFlightMode(sm, th);
                d = sm.descentSpeed;  v = sm.verticalSpinMag;  cn = sm.coneAngleDeg;
            else
                old = 'chaotic';            % the old grid's convention for a blow-up
            end
        catch
        end
        gLab{k} = lbl;  gOld{k} = old;  gParts{k} = prt;  gTilt(k) = tilt;  gRev(k) = rev;
        gDesc(k) = d;   gVs(k) = v;     gCone(k) = cn;
    end
    rs = @(v) reshape(v, nXg, cfg.gridN).';
    modeName  = string(rs(gLab));
    partsGrid = string(rs(gParts));
    [~, modeIdx]    = ismember(modeName, string(modeList));
    [~, modeIdxOld] = ismember(string(rs(gOld)), string(modeListOld));
    modeIdxOld(modeIdxOld == 0) = numel(modeListOld);
    gridMetrics = struct('descentSpeed',rs(gDesc), 'verticalSpinMag',rs(gVs), ...
                         'coneAngleDeg',rs(gCone), 'absSpanTiltDeg',rs(gTilt), ...
                         'revolving',rs(gRev));
    fprintf('  %d cells in %.0f s\n', total, toc(tGrid));

    gridTrack = summarizeModeGrid(modeName, gridMetrics.absSpanTiltDeg, gridMetrics.revolving);
    % q0/omega0 saved top-level too: the form planar/pickModeGridRuns.m reads, so
    % this file can be handed straight to the picker. modeList/modeCol are saved
    % so the picker draws the six-part classifier's names and colours.
    q0 = cfg.q0;   omega0 = cfg.omega0; %#ok<NASGU>
    save(fullfile(runDir,'mode_grid_3D.mat'), 'modeIdx','modeName','modeList','modeCol', ...
         'partsGrid','modeIdxOld','modeListOld','gridMetrics','gridTrack','chordFrac', ...
         'spanFrac','cfg','effSwitches','q0','omega0');

    summary(end+1) = "";
    summary(end+1) = "MODE GRID, our seed -- six-part classifier (rows = span offset desc, cols = chord offset)";
    summary(end+1) = "  codes: " + strjoin(compose('%s=%s', string(pal.codes(:)), ...
                                           string(modeList(:))), '  ');
    for iz = cfg.gridN:-1:1
        summary(end+1) = sprintf('  span %5.2f | %s', spanFrac(iz), pal.codes(modeIdx(iz,:)));
    end
    summary = [summary(:); gridTrack.lines];

    f = figure('Name','Mode grid 3D','Color','w','Visible','off','Position',[60 60 1500 650]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact');
    title(tl, sprintf('3D mode grid -- %s   (git %s)', cfg.label, gh), 'Interpreter','none');
    ax = nexttile(tl);
    imagesc(ax, spanFrac, chordFrac, modeIdx.');  set(ax,'YDir','normal');
    colormap(ax, modeCol);  caxis(ax, [0.5 numel(modeList)+0.5]);
    colorbar(ax, 'Ticks',1:numel(modeList), 'TickLabels',modeList);
    xlabel(ax,'spanwise nut offset (x span)');  ylabel(ax,'chordwise nut offset (x chord)');
    title(ax, sprintf('non-physical %.0f%%', 100*gridTrack.nonPhysicalFrac));
    ax = nexttile(tl);
    histogram(ax, 'BinEdges', gridTrack.tiltEdges, 'BinCounts', gridTrack.tiltHist);
    xline(ax, 30, 'k--', 'flat < 30');  xline(ax, 70, 'r--', 'end-on > 70');
    xline(ax, 11.4, 'b:', 'Hou AR 11.4');
    xlabel(ax, 'revolver |span tilt| from horizontal (deg)');  ylabel(ax, 'cells');
    title(ax, sprintf('%d revolvers, median %.1f deg', gridTrack.revolverCount, gridTrack.tiltMedian));
    drawnow;  exportgraphics(f, fullfile(runDir,'mode_grid_3D.png'), 'Resolution',150);
    close(f);
end

%% 8. Stage F -- PAPER COMPARISON (Hou et al. 2025)  [ON by default]
% Three labels per cell: the six-part classifier (of record, scored by
% summarizeModeGrid against the digitised Fig. 2a), the paper-vocabulary
% classifier (classifyPaperMode, with the 2026-10-03 sign fix), and the older
% classifyFlightMode.
paperTrack = [];  agreeFrac = NaN;
if cfg.runPaperGrid
    fprintf('\n--- F. paper comparison ---\n');
    [pcfg, pbsp, toNutPos, pinfo] = paperSeedConfig();
    for k = 1:numel(fn);  pcfg.(fn{k}) = sw.(fn{k});  end     % same overrides
    pcfg.tspan = cfg.tspan;  pcfg.odeRelTol = cfg.odeRelTol;  pcfg.odeAbsTol = cfg.odeAbsTol;
    pcfg.metricOpts = cfg.metricOpts;
    pmo = pcfg.metricOpts;
    pmo.refLength = pinfo.L;  pmo.sigma = pinfo.sigma;
    pmo.g = pcfg.g;  pmo.rhoFluid = pcfg.rhoFluid;
    pmo.spanLength = pinfo.L;                                 % their long axis = our span
    pmo.flipHysteresisTurns = cfg.seedModeThresholds.flipHysteresisTurns;
    pOde = odeset('RelTol',pcfg.odeRelTol,'AbsTol',pcfg.odeAbsTol);
    pq0 = pcfg.releaseQuat;  pw0 = pcfg.releaseOmega;
    pth = cfg.modeThresholds;  ths = cfg.seedModeThresholds;

    % --- the five published anchors ---------------------------------------
    anch = { 'AR ', 0.39, 0.17, -11.4,  4.2; ...
             'CST', 0.06, 0.04, -38.2,  2.3; ...
             'SST', 0.10, 0.10, -38.3, 13.0; ...
             'CH ', 0.18, 0.08,   NaN,  NaN; ...
             'FA ', 0.35, 0.08,   NaN,  NaN };
    summary(end+1) = "";
    summary(end+1) = "PAPER ANCHORS  (case | ours | paper-vocab | theta model | theta paper | V_d)";
    aLab = strings(size(anch,1),1);  aNew = aLab;
    for i = 1:size(anch,1)
        b = pbsp;  b.tSamples = 0;
        b.nutPos_t = toNutPos(anch{i,2}, anch{i,3});  b.nutMass_t = pcfg.nutMass;
        sp  = buildSeedParams(b, pcfg);
        rhs = seedRHS(sp);
        [t,x] = ode45(@(tt,xx) rhs(tt,xx,sp), pcfg.tspan, ...
                      [zeros(3,1); pq0; zeros(3,1); pw0], pOde);
        pm  = computeSeedModeMetrics(t, x, pmo);
        aLab(i) = string(classifyPaperMode(pm));
        aNew(i) = string(classifySeedMode(pm, ths));
        if isnan(anch{i,4}); tp = '(not reported)'; else; tp = sprintf('%+.1f+/-%.1f',anch{i,4},anch{i,5}); end
        summary(end+1) = sprintf('  %-4s | %-18s | %-4s | %+6.1f+/-%4.1f | %-14s | %.2f', ...
                                 anch{i,1}, aNew(i), aLab(i), pm.spanAxisTiltDeg, ...
                                 pm.spanAxisTiltStd, tp, pm.descentSpeed);
        fprintf('  anchor %-4s -> %-18s (%-4s)  theta %+.1f  V_d %.2f\n', anch{i,1}, aNew(i), ...
                aLab(i), pm.spanAxisTiltDeg, pm.descentSpeed);
    end

    % --- the phase map over their window ----------------------------------
    xaGrid = linspace(0.01, 0.41, cfg.paperNX);
    ybGrid = linspace(0.01, 0.20, cfg.paperNY);
    pTot = cfg.paperNX * cfg.paperNY;   pnX = cfg.paperNX;
    pLab = repmat({'failed'}, pTot, 1);  oLab = pLab;  nLab = pLab;
    pVd = nan(pTot,1);  pTh = nan(pTot,1);  pTilt = nan(pTot,1);  pRev = false(pTot,1);
    tP = tic;
    parfor k = 1:pTot
        iy = floor((k-1)/pnX)+1;   ix = mod(k-1,pnX)+1;
        b = pbsp;  b.tSamples = 0;
        b.nutPos_t = toNutPos(xaGrid(ix), ybGrid(iy));  b.nutMass_t = pcfg.nutMass;
        lab = 'failed';  our = 'failed';  nw = 'failed';
        vd = NaN;  thv = NaN;  tilt = NaN;  rev = false;
        try
            sp  = buildSeedParams(b, pcfg);
            rhs = seedRHS(sp);
            [t,x] = ode45(@(tt,xx) rhs(tt,xx,sp), pcfg.tspan, ...
                          [zeros(3,1); pq0; zeros(3,1); pw0], pOde);
            if all(isfinite(x(:)))
                pm  = computeSeedModeMetrics(t, x, pmo);
                lab = classifyPaperMode(pm);
                our = classifyFlightMode(pm, pth);
                [nw, p] = classifySeedMode(pm, ths);
                vd = pm.descentSpeed;  thv = pm.spanAxisTiltDeg;
                tilt = pm.absSpanTiltDeg;   rev = ~strcmp(p.revolution, 'none');
            end
        catch
        end
        pLab{k} = lab;  oLab{k} = our;  nLab{k} = nw;
        pVd(k) = vd;  pTh(k) = thv;  pTilt(k) = tilt;  pRev(k) = rev;
    end
    rp = @(v) reshape(v, pnX, cfg.paperNY).';
    paperMode = string(rp(pLab));   ourMode = string(rp(oLab));   newMode = string(rp(nLab));
    pVd = rp(pVd);   pTh = rp(pTh);   pTilt = rp(pTilt);   pRev = rp(pRev);
    refMode = strings(cfg.paperNY, pnX);
    for iy = 1:cfg.paperNY
        for ix = 1:pnX
            refMode(iy,ix) = string(paperModeReference(xaGrid(ix), ybGrid(iy)));
        end
    end
    % paper-vocabulary classifier agreement (kept for continuity with phase 5)
    fam = paperMode;  fam(fam=="CST" | fam=="SST") = "ST";
    valid = fam ~= "failed";
    agreeFrac = sum(fam(valid) == refMode(valid)) / max(nnz(valid),1);
    % six-part classifier: the tracked agreement, non-physical fraction, tilt
    paperTrack = summarizeModeGrid(newMode, pTilt, pRev, refMode);
    fprintf('  %d cells in %.0f s -- agreement %.0f%% (six-part), %.0f%% (paper-vocab)\n', ...
            pTot, toc(tP), 100*paperTrack.paperAgreeFrac, 100*agreeFrac);
    save(fullfile(runDir,'paper_mode_grid_3D.mat'), 'xaGrid','ybGrid','newMode','paperMode', ...
         'ourMode','refMode','pVd','pTh','pTilt','pRev','agreeFrac','paperTrack', ...
         'pcfg','pinfo','effSwitches');

    pal = seedModePalette();
    refOf = containers.Map({'AR','ST','CH','FA'}, {'A','T','X','F'});
    summary(end+1) = "";
    summary(end+1) = "PAPER GRID -- six-part classifier vs digitised Fig. 2a";
    summary = [summary(:); paperTrack.lines];
    summary(end+1) = sprintf('  (paper-vocabulary classifier, for continuity: %.0f%% agreement)', ...
                             100*agreeFrac);
    summary(end+1) = "  OURS (six-part codes)  THEIRS (A=AR T=ST X=CH F=FA)";
    for iy = cfg.paperNY:-1:1
        ro = blanks(pnX);  rr = blanks(pnX);
        for ix = 1:pnX
            j = find(strcmp(pal.names, newMode(iy,ix)), 1);
            ro(ix) = pal.codes(j);
            rr(ix) = refOf(char(refMode(iy,ix)));
        end
        summary(end+1) = sprintf('  y=%.3f | %-14s  | %s', ybGrid(iy), ro, rr);
    end
    summary(end+1) = "  six-part codes: " + strjoin(compose('%s=%s', string(pal.codes(:)), ...
                                                    string(pal.names(:))), '  ');
    summary(end+1) = sprintf('  descent %.2f-%.2f m/s (%.2f-%.2f x scale %.3f)', ...
        min(pVd(:)), max(pVd(:)), min(pVd(:))/pinfo.Vscale, max(pVd(:))/pinfo.Vscale, pinfo.Vscale);
end

%% 8b. Tracked numbers -- this run's file, plus one row in the running history
% tracking_3D.txt holds this run's numbers; tracking_history_3D.csv (in outPath)
% gets one row per run, so physics changes can be compared over time.
if ~isempty(gridTrack) || ~isempty(paperTrack)
    tr = strings(0,1);
    tr(end+1) = sprintf('3D test suite tracking  %s  git %s  label=%s', timestamp, gh, cfg.label);
    tr(end+1) = sprintf('physics: edgeDrag %d, rate %d, moment %d, LEV %d', ...
        effSwitches.enableEdgeDrag, effSwitches.enableAddedMassRate, ...
        effSwitches.enableAddedMassMoment, effSwitches.enableLEV);
    if ~isempty(gridTrack); tr = [tr(:); ""; "MODE GRID (our seed)"; gridTrack.lines]; end
    if ~isempty(paperTrack); tr = [tr(:); ""; "PAPER GRID (Hou et al. window)"; paperTrack.lines]; end
    fid = fopen(fullfile(runDir,'tracking_3D.txt'), 'w');
    fprintf(fid, '%s\n', tr);   fclose(fid);

    row = table(string(timestamp), string(gh), string(cfg.label), ...
        effSwitches.enableEdgeDrag, effSwitches.enableAddedMassRate, ...
        effSwitches.enableAddedMassMoment, effSwitches.enableLEV, ...
        'VariableNames', {'timestamp','git','label','edgeDrag','addedMassRate', ...
                          'addedMassMoment','LEV'});
    if ~isempty(gridTrack)
        row.grid_cells = gridTrack.nValid;            row.grid_nonPhysical = gridTrack.nonPhysicalFrac;
        row.grid_revolverFrac = gridTrack.revolverFrac; row.grid_tiltMedian = gridTrack.tiltMedian;
        row.grid_flatRevolverFrac = gridTrack.flatRevolverFrac;
    else
        row.grid_cells = NaN;  row.grid_nonPhysical = NaN;  row.grid_revolverFrac = NaN;
        row.grid_tiltMedian = NaN;  row.grid_flatRevolverFrac = NaN;
    end
    if ~isempty(paperTrack)
        row.paper_agreement = paperTrack.paperAgreeFrac;   row.paper_agreementPaperVocab = agreeFrac;
        row.paper_nonPhysical = paperTrack.nonPhysicalFrac; row.paper_tiltMedian = paperTrack.tiltMedian;
        row.paper_flatRevolverFrac = paperTrack.flatRevolverFrac;
    else
        row.paper_agreement = NaN;  row.paper_agreementPaperVocab = NaN;
        row.paper_nonPhysical = NaN;  row.paper_tiltMedian = NaN;  row.paper_flatRevolverFrac = NaN;
    end
    if isfile(cfg.trackingFile)
        writetable(row, cfg.trackingFile, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(row, cfg.trackingFile);
    end
    fprintf('tracked numbers -> %s (+1 row in %s)\n', fullfile(runDir,'tracking_3D.txt'), ...
            cfg.trackingFile);
end

%% 9. Manifest + summary
save(fullfile(runDir,'manifest_3D.mat'), 'cfg', 'sw', 'effSwitches', 'timestamp', 'gh', ...
     'gridTrack', 'paperTrack');
fid = fopen(fullfile(runDir,'summary_3D.txt'), 'w');
for k = 1:numel(summary);  fprintf(fid, '%s\n', summary(k));  end
fclose(fid);
fprintf('\n%s\n', strjoin(cellstr(summary), newline));
fprintf('\nWritten to %s\n', runDir);
