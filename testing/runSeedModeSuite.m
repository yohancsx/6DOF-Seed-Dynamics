%% Seed flight-MODE suite -- the seven named modes at their working inputs
% One section per mode, each a SINGLE drop using the exact nut position and
% initial condition from "Working Dynamics Inputs 7-24-26.txt", run with the
% working physics configuration. Each run prints a metric/classifier summary;
% the viewer at the bottom shows the 3D trajectory and Euler angles per mode so
% you can classify them by eye and debug toward the simplest model that still
% produces every mode.
%
% Body axes: x = chord, y = plate-normal, z = span.  c = chordLength, S = spanLength.
% Run Section 0 first (it resets `modeResults`), then the mode sections top to
% bottom, then the viewer.

%% 0. Configuration  -- working physics + base seed
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ are assumed already on the MATLAB path.

% --- Baseline seed geometry / material (working base seed) ----------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   % kg

% --- Environment ----------------------------------------------------------
cfg.rhoFluid = 1.225;   cfg.g = 9.81;

% --- Physics switches (the default "FULL-minus-geomVelocity" config) -------
cfg.enableSpanForce             = true;
cfg.enableSpanTorque            = true;
cfg.enableSpanGeomVelocity      = false;   % no measurable effect on any mode
cfg.enableSpanCOPMigration      = true;
cfg.enableSpanTorqueAttenuation = false;
cfg.enableTxDamping             = true;
cfg.aero = struct('C_span', 0.2, 'C_span_torque', 0.7);   % rest at defaults

% --- Simulation -----------------------------------------------------------
cfg.tspan = [0 10];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;

% --- Analysis (classifier is REFERENCE only; classify by eye from the plots)
cfg.metricOpts.windowStartFrac = 0.5;
cfg.metricOpts.convergeTol     = 0.20;
cfg.modeThresholds = defaultModeThresholds();

% --- Baseline base-seed-params --------------------------------------------
xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

% Shorthands + the pi/6-about-z release used by the fluttering modes.
c = cfg.chordLength;   S = cfg.spanLength;
qLevel = [1; 0; 0; 0];
qTilt  = axisAngleToQuat([0; 0; 1], pi/6);
noSpin = [0; 0; 0];

modeResults = [];   % reset the collection (Section 0 must run first)

%% 1. Spanwise-axis fluttering   nut center, pi/6 tilt about z
modeResults = [modeResults, ...
    runSingleMode('Spanwise-axis fluttering', [0; 0; 0], qTilt, noSpin, cfg, baseBsp)];

%% 2. Gliding                    nut +0.5c chordwise, from rest
modeResults = [modeResults, ...
    runSingleMode('Gliding', [0.5*c; 0; 0], qLevel, noSpin, cfg, baseBsp)];

%% 3. Diving                     nut +1.5c chordwise, from rest
modeResults = [modeResults, ...
    runSingleMode('Diving', [1.5*c; 0; 0], qLevel, noSpin, cfg, baseBsp)];

%% 4. Fluttering + spiral        nut +0.01S spanwise, pi/6 tilt about z
modeResults = [modeResults, ...
    runSingleMode('Fluttering + spiral', [0; 0; 0.01*S], qTilt, noSpin, cfg, baseBsp)];

%% 5. Fluttering + tight spiral  nut +S spanwise, pi/6 tilt about z
modeResults = [modeResults, ...
    runSingleMode('Fluttering + tight spiral', [0; 0; S], qTilt, noSpin, cfg, baseBsp)];

%% 6. Autorotation               nut +c chord & +1.2S span, from rest
modeResults = [modeResults, ...
    runSingleMode('Autorotation', [c; 0; 1.2*S], qLevel, noSpin, cfg, baseBsp)];

%% 7. Parachute                  nut -c below the plate (out-of-plane), from rest
modeResults = [modeResults, ...
    runSingleMode('Parachute', [0; -c; 0], qLevel, noSpin, cfg, baseBsp)];

%% 8. View each mode  -- trajectory + Euler angles
% Edit modesToView to a subset of indices (1..7) to view fewer at once.
modesToView = 1:numel(modeResults);
for i = modesToView
    r = modeResults(i);
    visualizeSeedTrajectory(r.t, r.x(:,1:3).', r.x(:,4:7).', ...
        struct('fig1Name', [r.name ' - traj'], 'fig2Name', [r.name ' - angles']));
end
