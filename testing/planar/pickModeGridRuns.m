%% Pick runs off a saved mode grid and animate them
% Interactive companion to runSeedModeGrid: open the phase map, CLICK the cells
% whose flight you want to watch, then batch-render a combined mode animation
% for each pick.
%
% MODEL-AGNOSTIC despite living in planar/: the re-integration dispatches through
% seedRHS on the saved cfg.shapeModel, so this works on a planar grid
% (Mode_Grid_results.mat) and equally on the 3D suite's mode_grid_3D.mat. Point
% resultsFile at either. It is the intended way to hand-check the classifier:
% click a cell, watch what the seed actually does, compare to the label.
%
% The grid only stores the mode LABEL per cell (the trajectories would be huge),
% so this script RE-INTEGRATES the dynamics at each selected nut position -- same
% cfg, release, and physics the grid used -- before handing (t, x) to
% animateModeTrajectory. Re-running is cheap (~2 s per cell).
%
% WORKFLOW (run as a live script, section by section):
%   Section 0 -> picks the grid: by default the NEWEST 3D-suite grid in
%                model_test_results/ (set resultsFile to open any other).
%   Section 1 -> a map window opens. LEFT-CLICK a cell to select it (a numbered
%                marker appears and its grid mode prints). Click as many as you
%                like. Press ANY KEY (or right-click) in the figure to finish.
%   Section 2 -> Run Section: one .mp4 per selected cell, written NEXT TO THE
%                GRID in <grid folder>/mode_grid_picks/, with a picks_log.csv
%                recording each pick (cell, label, parts, video).
%
% For grids labelled by the six-part classifier (the 3D suite's mode_grid_3D.mat),
% videos are BLIND (no mode colouring from the older classifier, which would
% contradict the label) and titled with the grid's six-part label and parts.
%
% Body axes: x = chord, y = normal, z = span.

%% 0. Configuration  -- EDIT here if needed
root = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% Grid to pick from. "" = the NEWEST model_test_results/*_3D_*/mode_grid_3D.mat
% (run folders start with a timestamp, so the last by name is the newest). Set a
% path to open any other grid, including a planar Mode_Grid_results.mat.
resultsFile = "";
if resultsFile == ""
    g = dir(fullfile(root, 'model_test_results', '*_3D_*', 'mode_grid_3D.mat'));
    if isempty(g)
        error('pickModeGridRuns:noGrid', 'No mode_grid_3D.mat under %s.', ...
              fullfile(root, 'model_test_results'));
    end
    [~, o] = sort({g.folder});   g = g(o);
    resultsFile = string(fullfile(g(end).folder, g(end).name));
end
fprintf('Grid: %s\n', resultsFile);

%% 1. Open the phase map and click the runs to animate
D = load(resultsFile);
cfg = D.cfg;   chordFrac = D.chordFrac;   spanFrac = D.spanFrac;   modeIdx = D.modeIdx;
[Nz, Nx] = size(modeIdx);          % Nz span rows, Nx chord cols
% A grid saved with its own vocabulary (the 3D suite saves the six-part
% classifier's modeList/modeCol) is drawn with it; older grids use the original
% classifyFlightMode palette.
if isfield(D, 'modeList') && isfield(D, 'modeCol')
    modeList = D.modeList;   modeCol = D.modeCol;
else
    [modeList, modeCol] = seedModeColors();
end

% Release condition + base seed -- reconstruct EXACTLY what the grid ran. Use the
% saved release if present (newer grids), else the grid's default pi/6-about-z.
if isfield(D, 'q0')     && ~isempty(D.q0);     q0 = D.q0;         else; q0 = axisAngleToQuat([0;0;1], pi/6); end
if isfield(D, 'omega0') && ~isempty(D.omega0); omega0 = D.omega0; else; omega0 = [0;0;0]; end
xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

% Map axes: span horizontal, chord vertical, both in chord units (as in the grid).
sc = cfg.spanLength / cfg.chordLength;
xSpan = spanFrac * sc;   yChord = chordFrac;
MI  = modeIdx.';                                   % Nx (chord rows) x Nz (span cols)
img = cat(3, reshape(modeCol(MI,1),Nx,Nz), reshape(modeCol(MI,2),Nx,Nz), reshape(modeCol(MI,3),Nx,Nz));

figPick = figure('Name', 'Pick mode-grid runs', 'Color', 'w');
image(xSpan, yChord, img);   set(gca, 'YDir', 'normal');   axis image;   hold on;
hs = sc/2;   hc = 0.5;
plot([-hs hs hs -hs -hs], [-hc -hc hc hc -hc], 'k--', 'LineWidth', 1.2);   % seed body
% one legend entry per mode present
present = unique(modeIdx(:)).';
for j = present
    plot(nan, nan, 's', 'MarkerFaceColor', modeCol(j,:), 'MarkerEdgeColor', 'k', ...
         'MarkerSize', 10, 'LineStyle', 'none', 'DisplayName', modeList{j});
end
legend('Location', 'eastoutside');
xlabel('Spanwise nut offset  (\times chord)');   ylabel('Chordwise nut offset  (\times chord)');
title('LEFT-CLICK cells to select  --  ANY KEY / right-click to finish');
xlim([min(xSpan) max(xSpan)]);   ylim([min(yChord) max(yChord)]);

% --- Click loop: left-click adds a (snapped) cell, any key / right-click ends --
sel = zeros(0, 2);                 % rows of [iz ix]
fprintf('Click runs on the map (left-click add, any key / right-click to finish)...\n');
while true
    [xc, yc, btn] = ginput(1);
    if isempty(btn) || btn ~= 1    % Enter (empty), a keypress, or a non-left click
        break;
    end
    [~, iz] = min(abs(xSpan  - xc));      % nearest span cell
    [~, ix] = min(abs(yChord - yc));      % nearest chord cell
    if ~isempty(sel) && any(sel(:,1)==iz & sel(:,2)==ix)
        continue;                          % already picked -> ignore
    end
    sel(end+1, :) = [iz ix];                                   %#ok<SAGROW>
    k = size(sel, 1);
    plot(xSpan(iz), yChord(ix), 'ko', 'MarkerSize', 11, 'LineWidth', 1.5, ...
         'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    text(xSpan(iz), yChord(ix), sprintf(' %d', k), 'FontWeight', 'bold', ...
         'VerticalAlignment', 'middle');
    fprintf('  [%d] chord=%.2f  span=%.2f  ->  %s\n', ...
            k, chordFrac(ix), spanFrac(iz), modeList{modeIdx(iz,ix)});
end
fprintf('Selected %d run(s). Now run Section 2.\n', size(sel,1));

%% 2. Generate animations for the selected runs
% Written beside the grid file, so a run folder holds its grid AND the picks
% made from it.
outputFolder = fullfile(fileparts(resultsFile), "mode_grid_picks");
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

% Animation options (forwarded to animateModeTrajectory). playbackSpeed<1 = slow-mo.
animOpts = struct('fps', 30, 'playbackSpeed', 1, 'showSeedVels', false);

if isempty(sel)
    warning('No runs selected -- run Section 1 and click some cells first.');
else
    fprintf('Rendering %d animation(s) to %s\n', size(sel,1), outputFolder);
    for s = 1:size(sel, 1)
        fprintf('(%d/%d) ', s, size(sel,1));
        animateGridRun(sel(s,1), sel(s,2), D, baseBsp, q0, omega0, outputFolder, animOpts);
    end
    fprintf('Done.\n');
end


% =========================================================================
% LOCAL: re-integrate one grid cell and render its combined mode animation
% =========================================================================
function animateGridRun(iz, ix, D, baseBsp, q0, omega0, outputFolder, animOpts)
    cfg = D.cfg;   c = cfg.chordLength;   S = cfg.spanLength;
    if isfield(D, 'modeList'); modeList = D.modeList; else; [modeList, ~] = seedModeColors(); end

    % Rebuild the seed at this nut position and integrate (same as the grid).
    bsp = baseBsp;
    bsp.tSamples  = cfg.tSamples;
    bsp.nutPos_t  = [D.chordFrac(ix)*c; 0; D.spanFrac(iz)*S];
    bsp.nutMass_t = cfg.nutMass;
    sp = buildSeedParams(bsp, cfg);
    x0 = [zeros(3,1); q0(:); zeros(3,1); omega0(:)];
    odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);
    % seedRHS, not a hardcoded seed6DOFODE: the saved cfg carries its shapeModel,
    % so a grid produced under shape3d re-integrates under the 3D RHS and a planar
    % grid under the planar one. Hardcoding the planar RHS here would have silently
    % animated DIFFERENT physics from the one that labelled the cell.
    rhs = seedRHS(sp);
    [t, x] = ode45(@(tt,xx) rhs(tt, xx, sp), cfg.tspan, x0, odeOpts);

    lbl = modeList{D.modeIdx(iz,ix)};
    sixPart = isfield(D, 'partsGrid');            % labelled by classifySeedMode
    parts = '';
    if sixPart; parts = char(D.partsGrid(iz,ix)); end
    tag = sprintf('chord%03.0f_span%03.0f_%s', 100*D.chordFrac(ix), 100*D.spanFrac(iz), lbl);
    vfile = fullfile(outputFolder, [matlab.lang.makeValidName(tag) '.mp4']);
    ttl = sprintf('nut chord=%.2fc  span=%.2fS   (grid mode: %s)', ...
                  D.chordFrac(ix), D.spanFrac(iz), lbl);
    if sixPart; ttl = sprintf('%s  [%s]', ttl, parts); end

    o = animOpts;   o.videoFile = char(vfile);   o.title = ttl;   o.seedParams = sp;
    if sixPart
        % no colouring by the older classifier's timeline -- it would contradict
        % the six-part label in the title; the ground track shows any revolution
        o.blind = true;   o.topPanel = 'groundTrack';
    end
    try
        animateModeTrajectory(t, x, o);
        close(gcf);
        fprintf('    wrote %s\n', vfile);
    catch ME
        if ~isempty(findobj('Type','figure')); close(gcf); end
        fprintf(2, '    FAILED %s: %s\n', tag, ME.message);
        return
    end

    % one row per pick, so the folder records what was looked at and why
    row = table(string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                D.chordFrac(ix), D.spanFrac(iz), string(lbl), string(parts), ...
                string([matlab.lang.makeValidName(tag) '.mp4']), ...
                'VariableNames', {'time','chordFrac','spanFrac','gridLabel','parts','video'});
    logFile = fullfile(outputFolder, 'picks_log.csv');
    if isfile(logFile)
        writetable(row, logFile, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(row, logFile);
    end
end
