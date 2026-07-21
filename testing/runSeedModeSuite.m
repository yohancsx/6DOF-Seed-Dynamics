%% Seed flight-MODE suite -- elicit & classify biological descent modes
% Five sections, one per target mode. Each moves the nut along the relevant
% axis over a tunable range, drops the seed, and uses computeTrajectoryMetrics
% + classifyFlightMode to name the resulting flight mode -- so you can (a) tune
% the sweep until the mode appears clearly, and (b) check that the classifier
% labels kinematics you already recognise.
%
% Adjust the per-section `sweep.*` settings (range, increments, axis) and the
% classifier thresholds in Section 0 as you go. Run section-by-section.
%
% Body axes: x=chord, y=normal/vertical, z=span. Sweep axes:
%   chord    -> glide / dive        span  -> spiral tumbling / tight spiral
%   diag     -> autorotation        vertical (out-of-plane Y) -> parachuting
%
% Free of local functions, so it can be saved as a live script (.mlx).

%% 0. Configuration & baseline seed  -- EDIT HERE
% --- Paths (edit for your machine) ----------------------------------------
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ are assumed already on the MATLAB path.

% --- Baseline seed geometry / material ------------------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;    cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;    % nut mass (kg); its POSITION is what each section sweeps

% --- Environment ----------------------------------------------------------
cfg.rhoFluid = 1.225;   cfg.g = 9.81;   cfg.enableSpanForce = false;

% --- Simulation (modes need time to develop) ------------------------------
cfg.tspan = [0 8];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;

% --- Analysis -------------------------------------------------------------
cfg.metricOpts.windowStartFrac = 0.5;    % use the latter 50% for steady metrics
cfg.metricOpts.convergeTol     = 0.20;
cfg.modeThresholds = defaultModeThresholds();   % <-- tune these fields to taste
cfg.plotVisible    = 'on';               % show overlay figures live

% --- Baseline base-seed-params (symmetric rectangular reference) ----------
xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baselineBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baselineBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baselineBsp.seedThickness = cfg.thickness;
baselineBsp.numStrips     = cfg.numStrips;
baselineBsp.tSamples      = cfg.tSamples;
baselineBsp.nutMass_t     = cfg.nutMass * ones(size(cfg.tSamples));
baselineBsp.nutPos_t      = repmat([0;0;0], 1, numel(cfg.tSamples));

baselineSeedParams = buildSeedParams(baselineBsp, cfg);

%% 1. Gliding / diving  (nut along the CHORD)
% Increasing chordwise offset should pass flutter -> glide -> dive. Widen the
% range or add increments until you clearly see a glide and then a dive.
sweep = struct();
sweep.name         = 'glide_dive';
sweep.axis         = 'chord';
sweep.fracRange    = [0.0 2.0];   % 0 .. 2x half-chord (off-plate at >1)
sweep.nSweep       = 6;
sweep.expectedMode = {'fluttering', 'gliding', 'diving'};
res_glideDive = runModeSweep(sweep, cfg, baselineBsp, baselineSeedParams);

%% 2. Spiral tumbling  (nut along the SPAN, small-to-intermediate)
% A modest spanwise offset should tumble about the spanwise axis while circling
% the vertical axis -- a wide helix.
sweep = struct();
sweep.name         = 'spiral_tumbling';
sweep.axis         = 'span';
sweep.fracRange    = [0.2 1.0];
sweep.nSweep       = 6;
sweep.expectedMode = {'spiralTumbling', 'tumbling'};
res_spiralTumbling = runModeSweep(sweep, cfg, baselineBsp, baselineSeedParams);

%% 3. Tight spiral  (nut along the SPAN, further out)
% Larger spanwise offset should tighten the helix; note it may transition
% toward autorotation as the seed becomes strongly tip-heavy.
sweep = struct();
sweep.name         = 'tight_spiral';
sweep.axis         = 'span';
sweep.fracRange    = [1.0 2.0];
sweep.nSweep       = 6;
sweep.expectedMode = {'tightSpiral', 'autorotation'};
res_tightSpiral = runModeSweep(sweep, cfg, baselineBsp, baselineSeedParams);

%% 4. Autorotation  (nut along the DIAGONAL, off-body; small spin nudge)
% Strong offset (nut outside the body) plus a tiny yaw nudge to break the
% unstable from-rest symmetry -> steady spin about the vertical axis.
sweep = struct();
sweep.name         = 'autorotation';
sweep.axis         = 'diag';
sweep.fracRange    = [1.0 2.5];
sweep.nSweep       = 6;
sweep.expectedMode = {'autorotation', 'tightSpiral'};
sweep.omega0       = [0; 2; 0];   % small initial spin about the plate normal (rad/s)
res_autorotation = runModeSweep(sweep, cfg, baselineBsp, baselineSeedParams);

%% 5. Parachuting  (nut BELOW the plate, out-of-plane along Y)
% CoM below the plate is pendulum-stable broadside -> slow, near-vertical,
% non-spinning descent. First exercise of a nonzero-Y nut -- sanity-check it.
sweep = struct();
sweep.name         = 'parachuting';
sweep.axis         = 'vertical';
sweep.fracRange    = [1.0 5.0];   % nut 1..5 half-chords below the plate
sweep.nSweep       = 5;
sweep.expectedMode = {'parachuting'};
res_parachuting = runModeSweep(sweep, cfg, baselineBsp, baselineSeedParams);
