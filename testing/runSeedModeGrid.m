%% Flight-mode phase map -- a Fig-2-style map of mode vs CoM (nut) position
% Sweeps a grid of nut positions in the (chordwise, spanwise) plane, drops the
% seed at each, classifies the resulting flight mode, and renders a coloured map
% of where each mode appears -- the simulation analogue of Fig. 2 of "Aerodynamic
% significance of mass distribution on diverse samara descent behaviors".
%
% The spanwise axis uses a NON-UNIFORM grid, dense near zero, because the
% tumbling+spiral band sits at very small spanwise offsets. Edit chordFrac /
% spanFrac in Section 0 to zoom or refine. Uses the default physics config.
%
% Body axes: x = chord, y = normal, z = span. Grid axes are nut offset as a
% fraction of chordLength (x) and spanLength (z).

%% 0. Configuration  -- EDIT HERE
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ assumed on the MATLAB path.

cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   cfg.rhoFluid = 1.225;   cfg.g = 9.81;
cfg.tspan = [0 12];   cfg.odeRelTol = 1e-5;   cfg.odeAbsTol = 1e-7;   % longer so modes settle; loose tol keeps the grid quick
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();
% Physics: leave switches/aero UNSET so the code defaults (FULL-minus-geomVelocity) apply.

% --- The grid (nut offset as fractions of chord / span) -------------------
% Modest default so it runs in a few minutes; raise the counts for a finer map.
chordFrac = linspace(0, 1.5, 40);                       % chordwise: uniform 0..2 c
% spanFrac  = [0 0.01 0.02 0.05 0.08 0.10 0.16 0.24 ...             % spanwise: dense near 0...
%              0.35 0.50 0.70 0.90 1.10 1.5 2];                 % ...then out to 1.1 S
spanFrac = linspace(0, 1.5, 40);
% --- Release condition (pi/6 tilt about z, as the spiral modes need) ------
q0     = axisAngleToQuat([0; 0; 1], pi/6);
omega0 = [0; 0; 0];

% --- Plot options ---------------------------------------------------------
showMarkers = true;   % false -> just the interpolated colour field (no per-cell markers)

% --- Baseline base-seed-params --------------------------------------------
xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

c = cfg.chordLength;   S = cfg.spanLength;
odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);

%% 1. Run the grid  (build seed -> integrate -> classify, per cell)
Nx = numel(chordFrac);   Nz = numel(spanFrac);
modeGrid = cell(Nz, Nx);                 % mode label per cell
fprintf('Mode grid: %d chord x %d span = %d runs\n', Nx, Nz, Nx*Nz);
tStart = tic;
for iz = 1:Nz
    for ix = 1:Nx
        nutPos = [chordFrac(ix)*c; 0; spanFrac(iz)*S];
        bsp = baseBsp;
        bsp.tSamples  = cfg.tSamples;
        bsp.nutPos_t  = nutPos;
        bsp.nutMass_t = cfg.nutMass;
        sp = buildSeedParams(bsp, cfg);
        x0 = [zeros(3,1); q0; zeros(3,1); omega0];
        try
            [t, x] = ode45(@(tt,xx) seed6DOFODE(tt, xx, sp), cfg.tspan, x0, odeOpts);
            if any(~isfinite(x(:)))
                modeGrid{iz, ix} = 'chaotic';       % diverged
            else
                m = computeTrajectoryMetrics(t, x, cfg.metricOpts);
                modeGrid{iz, ix} = classifyFlightMode(m, cfg.modeThresholds);
            end
        catch
            modeGrid{iz, ix} = 'failed';            % integrator error
        end
    end
    fprintf('  span row %2d/%2d (frac %.2f) done  [%.0f s]\n', iz, Nz, spanFrac(iz), toc(tStart));
end

%% 2. Mode legend (label -> integer index -> colour) + text summary
modeList = {'gliding','diving','parachuting','fluttering','spiral', ...
            'tightSpiral','autorotation','chaotic','undetermined','failed'};
modeCode = 'gdpfstacux';                 % single-letter codes for the ASCII summary
modeCol  = [0.20 0.70 0.20;   % gliding      green
            0.70 0.10 0.10;   % diving       dark red
            0.20 0.80 0.80;   % parachuting  cyan
            0.95 0.60 0.10;   % fluttering   orange
            0.20 0.40 0.90;   % spiral       blue
            0.55 0.20 0.70;   % tightSpiral  purple
            0.92 0.20 0.60;   % autorotation magenta
            0.30 0.30 0.30;   % chaotic      dark grey
            0.75 0.75 0.75;   % undetermined light grey
            0.00 0.00 0.00];  % failed       black

modeIdx = zeros(Nz, Nx);
for iz = 1:Nz
    for ix = 1:Nx
        j = find(strcmp(modeGrid{iz, ix}, modeList), 1);
        if isempty(j); j = numel(modeList); end   % unknown -> 'failed'
        modeIdx(iz, ix) = j;
    end
end

% Console ASCII map (span rows top = large offset, so it reads like the plot).
fprintf('\nMode map (rows=span offset, cols=chord offset). Codes: %s\n', ...
        strjoin(cellfun(@(s,c)[c '=' s], modeList, num2cell(modeCode), 'uni', 0), '  '));
for iz = Nz:-1:1
    fprintf('span %5.2f | %s\n', spanFrac(iz), modeCode(modeIdx(iz, :)));
end
[uc, ~, ui] = unique(modeIdx(:));
fprintf('\ncounts: ');
for k = 1:numel(uc); fprintf('%s=%d  ', modeList{uc(k)}, sum(ui==k)); end
fprintf('\n');

%% 3. Coloured phase map + save
% Axes: SPAN offset horizontal, CHORD offset vertical, both expressed in units
% of CHORD (so span uses the span/chord ratio). Equal data aspect makes the plot
% read at the seed's true proportions (span is span/chord times longer than the
% chord), and a dashed rectangle outlines the actual seed body.
sc = cfg.spanLength / cfg.chordLength;   % span-per-chord (physical aspect factor)
xSpan  = spanFrac  * sc;                 % horizontal: span offset in chord units
yChord = chordFrac;                      % vertical:   chord offset in chord units

% Per-cell RGB, transposed to (chord rows) x (span cols) for the swapped axes.
MI = modeIdx.';                          % Nx (chord) x Nz (span)
cr = modeCol(:, 1);   cg = modeCol(:, 2);   cb = modeCol(:, 3);
R = cr(MI);   G = cg(MI);   B = cb(MI);

xf = linspace(min(xSpan),  max(xSpan),  300);     % fine uniform display grid
yf = linspace(min(yChord), max(yChord), 300);
[XF, YF] = meshgrid(xf, yf);
rgbImg = cat(3, interp2(xSpan, yChord, R, XF, YF, 'linear'), ...
                interp2(xSpan, yChord, G, XF, YF, 'linear'), ...
                interp2(xSpan, yChord, B, XF, YF, 'linear'));
rgbImg = min(max(rgbImg, 0), 1);                  % clamp to valid RGB

figure('Name', 'Flight-mode phase map', 'Color', 'w');
image(xf, yf, rgbImg);           % interpolated colour background
set(gca, 'YDir', 'normal');      % chord offset increases upward
hold on;

% Legend + optional markers: one entry per mode present. Markers carry a black
% edge so they stay visible on the coloured field.
[XX, YY] = meshgrid(xSpan, yChord);   % Nx x Nz, matching MI
present = unique(MI(:)).';
for j = present
    if showMarkers
        sel = (MI == j);
        scatter(XX(sel), YY(sel), 55, 'filled', 's', ...
                'MarkerFaceColor', modeCol(j, :), 'MarkerEdgeColor', 'k', ...
                'LineWidth', 0.4, 'DisplayName', modeList{j});
    else
        plot(nan, nan, 's', 'MarkerFaceColor', modeCol(j, :), 'MarkerEdgeColor', 'k', ...
             'MarkerSize', 10, 'LineStyle', 'none', 'DisplayName', modeList{j});
    end
end

% Seed-body outline (dashed black), drawn LAST so it sits on top. The seed is
% centred at the nut-offset origin and extends +/- half a span / half a chord;
% in these chord-units that is +/- sc/2 horizontally and +/- 0.5 vertically.
% Only the positive quadrant is in view, so its top and right edges are shown.
hs = sc / 2;   hc = 0.5;
plot([-hs hs hs -hs -hs], [-hc -hc hc hc -hc], 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'seed body');

axis equal;                       % true seed proportions
box on;
xlabel('Spanwise nut offset  (\times chord)');
ylabel('Chordwise nut offset  (\times chord)');
title('Flight mode vs. CoM (nut) position');
legend('Location', 'eastoutside');
xlim([min(xSpan),  max(xSpan)]);
ylim([min(yChord), max(yChord)]);

outDir = char(fileparts(char(helpersFolder)));
savefig(gcf, fullfile(outDir, 'Mode_Grid.fig'));
exportgraphics(gcf, fullfile(outDir, 'Mode_Grid.png'), 'Resolution', 150);
save(fullfile(outDir, 'Mode_Grid_results.mat'), 'chordFrac', 'spanFrac', 'modeGrid', 'modeIdx', 'cfg');
fprintf('\nSaved Mode_Grid.fig / .png / _results.mat to %s\n', outDir);
