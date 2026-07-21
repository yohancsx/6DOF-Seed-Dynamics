%% Seed dynamics test suite -- one-click sweep runner
% Runs a battery of seed-dynamics simulations (nut mass/position sweeps, initial
% pitch/roll tilts, initial yaw spin, one-sided lift/drag asymmetry, and a
% strip-count convergence check), saves a 3D-trajectory .fig + .png and a
% loadable .mat per case, and a colour-graded overlay figure per group.
%
% This is NOT validation against ground truth -- it is a consistent, repeatable
% picture of what the model does, so that after changing the model you can re-run
% and eyeball what shifted. Each run is written to a timestamped subfolder.
%
% Cases default to "dropped from rest, level, no spin" unless the case specifies
% otherwise. Body axes: x=chord, y=normal/vertical, z=span. See seedTestCases.m.
%
% To re-plot a single case later (e.g. its Euler angles), load its .mat and call:
%   load('.../yawSpin_3.0rads.mat');           % gives struct `result`
%   visualizeSeedTrajectory(result.t, result.x(:,1:3).', result.x(:,4:7).');
%
% This script is section-organised and free of local functions, so it can be
% saved as a MATLAB live script (.mlx) unchanged.

%% 1. Configuration  -- EDIT HERE
% --- Paths (edit for your machine) ----------------------------------------
% Where all run outputs are written:
outPath       = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\test_suite";
% Folder holding this suite's helper functions (seedTestCases, runOneSeedCase, ...):
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: the model code (physics/ and visualization/) is assumed to be on the
% MATLAB path already (add it once per session if not).

% --- Output location ------------------------------------------------------
cfg.outputRoot = outPath;

% --- Baseline seed geometry / material ------------------------------------
cfg.spanLength  = 0.050;   % spanwise length, body z (m)
cfg.chordLength = 0.015;   % chord length,    body x (m)
cfg.thickness   = 0.002;   % plate thickness        (m)
cfg.bulkDensity = 65;    % bulk density           (kg/m^3)
cfg.numStrips   = 10;      % baseline strip count
cfg.tSamples    = 0;       % mass-props time samples (constant mass -> scalar)
cfg.nutMass     = 75e-6;   % baseline nut mass (kg)
cfg.nutPos      = [0;0;0]; % baseline nut position, body [x;y;z] (m)

% --- Environment ----------------------------------------------------------
cfg.rhoFluid        = 1.225;   % air density (kg/m^3)
cfg.g               = 9.81;    % gravity (m/s^2)
cfg.enableSpanForce = true;   % whole-seed spanwise force (on by default)
% cfg.aero = struct(...);      % OPTIONAL computeAeroCoeffs overrides (else defaults)

% --- Simulation -----------------------------------------------------------
cfg.tspan     = [0 5];    % drop duration (s) -- some dynamics take ~5 s to appear
cfg.odeRelTol = 1e-6;
cfg.odeAbsTol = 1e-8;

% --- Sweep settings (all editable) ----------------------------------------
cfg.nIncr         = 5;          % increments per numeric sweep
cfg.nutMassFrac   = 0.20;       % nut mass +/- 20%
cfg.nutPosMaxFrac = 1.20;       % nut position 0 .. 120% of half-dimension (>1 = off-plate, allowed)
cfg.tiltMaxDeg    = 45;         % pitch / roll drop from -45 .. +45 deg
cfg.yawSpinMin    = 1;          % initial yaw spin sweep (rad/s)
cfg.yawSpinMax    = 5;
cfg.asymFactor    = 0.5;        % one-sided lift/drag multiplier
cfg.stripCounts   = [1 5 10 20];% strip-convergence counts

% --- Which groups to run (toggle to run a subset) -------------------------
cfg.groups.nutMass   = true;
cfg.groups.nutChord  = true;
cfg.groups.nutSpan   = true;
cfg.groups.nutDiag   = true;
cfg.groups.pitch     = true;
cfg.groups.roll      = true;
cfg.groups.yawSpin   = true;
cfg.groups.asymmetry = true;
cfg.groups.stripConv = true;

% --- Output options -------------------------------------------------------
cfg.savePng     = true;    % also save .png next to each .fig
cfg.saveOverlay = true;    % save per-group overlay figures

%% 2. Build the case list
[cases, baselineBsp] = seedTestCases(cfg);
fprintf('Test suite: %d cases across enabled groups.\n', numel(cases));

%% 3. Prepare the run output folder
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
runDir    = fullfile(cfg.outputRoot, timestamp);
if ~exist(runDir, 'dir'); mkdir(runDir); end
fprintf('Output: %s\n', runDir);

%% 4. Build the baseline seed once (reused by initial-condition-only cases)
baselineSeedParams = buildSeedParams(baselineBsp, cfg);

%% 5. Run every case: integrate, plot, save
results = [];
for k = 1:numel(cases)
    cs = cases(k);
    fprintf('[%2d/%2d] %-10s %-22s ... ', k, numel(cases), cs.group, cs.label);

    r = runOneSeedCase(cs, baselineSeedParams, baselineBsp, cfg);
    fprintf('%s\n', r.status);

    % Per-group subfolder.
    caseDir = fullfile(runDir, r.group);
    if ~exist(caseDir, 'dir'); mkdir(caseDir); end

    % Loadable per-case result (t, x, x0, seedParams, status, ...).
    result = r; %#ok<NASGU>
    save(fullfile(caseDir, [r.label '.mat']), 'result');

    % Per-case 3D trajectory figure (whenever there is a trajectory to draw).
    if ~isempty(r.x)
        fig = plotSeedCaseTrajectory(r, 'off');
        savefig(fig, fullfile(caseDir, [r.label '.fig']));
        if cfg.savePng
            exportgraphics(fig, fullfile(caseDir, [r.label '.png']), 'Resolution', 150);
        end
        close(fig);
    end

    results = [results; r]; %#ok<AGROW>
end

%% 6. Per-group overlay figures
if cfg.saveOverlay
    groupNames = unique({results.group}, 'stable');
    for gi = 1:numel(groupNames)
        gName = groupNames{gi};
        gRes  = results(strcmp({results.group}, gName));
        fig   = plotGroupOverlay(gRes, gName, 'off');
        savefig(fig, fullfile(runDir, gName, '_overlay.fig'));
        if cfg.savePng
            exportgraphics(fig, fullfile(runDir, gName, '_overlay.png'), 'Resolution', 150);
        end
        close(fig);
    end
end

%% 7. Manifest + summary
% Compact per-case summary (no bulky trajectories) for the manifest + table.
caseSummary = struct('group',     {results.group}, ...
                     'label',     {results.label}, ...
                     'status',    {results.status}, ...
                     'sweepVal',  {results.sweepVal}, ...
                     'descended', {results.descended});
save(fullfile(runDir, 'manifest.mat'), 'cfg', 'timestamp', 'caseSummary');

% Human-readable summary table (also echoed to the console).
summaryPath = fullfile(runDir, 'summary.txt');
fid = fopen(summaryPath, 'w');
header = sprintf('Seed test suite  %s\n\n%-11s %-22s %-10s %-9s\n', ...
                 timestamp, 'group', 'label', 'status', 'descended');
fprintf('%s', header);
fprintf(fid, '%s', header);
for k = 1:numel(results)
    line = sprintf('%-11s %-22s %-10s %-9d\n', results(k).group, results(k).label, ...
                   results(k).status, results(k).descended);
    fprintf('%s', line);
    fprintf(fid, '%s', line);
end
fclose(fid);

nBad = sum(~strcmp({results.status}, 'OK'));
fprintf('\nDone. %d cases, %d not OK. Results in:\n  %s\n', numel(results), nBad, runDir);
