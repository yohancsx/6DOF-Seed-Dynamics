%% Seed Dynamics ODE Test -- 3D MODEL
% The shape3d counterpart of Seed_Dynamics_ODE_Test.m: a single-drop harness for
% quick, hand-driven experiments. Edit section 1, run the file, look at the
% trajectory and the video. Nothing here sweeps, scores or saves a run folder --
% that is runSeedTestSuite3D.m's job. This is the scratchpad.
%
% WHAT IT ADDS OVER THE 2D VERSION
%   - builds through setupSeedShape3D (via buildSeedParams) and integrates with
%     seed6DOFODE3D, dispatched by seedRHS on cfg.shapeModel;
%   - TWIST and CURVATURE profiles, the two shape mechanisms that exist only in
%     this model, both off by default so the seed starts flat;
%   - the shape3d physics switches (edge drag, the two Kirchhoff added-mass
%     couplings, LEV) in one editable block, applied through buildSeedParams so a
%     switch this model does not honour WARNS instead of sitting there inert;
%   - a metrics + classifier line printed after the run, so a case can be checked
%     against the label the sweeps would have given it.
%
% Body axes: x = chord, y = plate normal, z = span. World Y is up.
% Sections: 1. seed + environment  2. shape sanity check  3. integrate
%           4. trajectory + metrics  5. animate and save
% Section-organised and free of local functions, so it saves as an .mlx cleanly.

%% 1. Setup seed shape, mass, environment and physics switches  -- EDIT HERE
root = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% --- Geometry & material --------------------------------------------------
spanLength   = 0.050;   % spanwise length, body z (m)
chordLength  = 0.015;   % chord length,    body x (m)
thickness    = 0.002;   % plate thickness        (m)
bulkDensity  = 65;      % bulk material density  (kg/m^3)
numStrips    = 20;      % spanwise aerodynamic strips
arealDensity = bulkDensity * thickness;
hs = spanLength / 2;    hc = chordLength / 2;   % half-span, half-chord

% --- Nut (the lumped seed mass) -------------------------------------------
nutMass = 75e-6;          % kg
nutPos  = [0; 0; 0];      % body [x;y;z] (m). e.g. [0.5*hc; 0; 0.5*hs] to offset it
tSamples = 0;             % single sample -> constant mass. A vector here (with a
                          % 3xN nutPos_t) drives the time-varying mass path.

% --- SHAPE: twist and curvature (shape3d only; [] = flat) -----------------
% twist(z)     : geometric pitch of each strip about its OWN span axis (rad),
%                as a function of spanwise station z (m). Strip centres stay at
%                y = 0, so a twisted seed is still planar in the mass sense.
% curvature(s) : dihedral angle (rad) vs span arc length s (m). This one lifts
%                the strips OUT of the plane, so it moves the CoM as well.
% Both take and return vectors. Comment out for a flat seed.
tipTwistDeg = 0;    % 0 -> no twist
tipCurveDeg = 0;    % 0 -> no curvature
twistFn     = [];   curvatureFn = [];
if tipTwistDeg ~= 0
    twistFn     = @(z) (deg2rad(tipTwistDeg)/hs) .* z;          % anti-symmetric
end
if tipCurveDeg ~= 0
    curvatureFn = @(s) (deg2rad(tipCurveDeg)/hs) .* s;          % symmetric bowl
    % one side only:  @(s) (s > 0) .* (deg2rad(tipCurveDeg)/hs) .* s
end

% --- Environment ----------------------------------------------------------
cfg.rhoFluid   = 1.225;   % air density (kg/m^3)
cfg.g          = 9.81;    % gravity (m/s^2)
cfg.shapeModel = 'shape3d';
% cfg.aero = struct(...);   % OPTIONAL computeAeroCoeffs overrides (else defaults)

% --- PHYSICS SWITCHES -- comment a line out to keep the builder default ---
% Defaults (setupSeedShape3D): edge drag ON, added-mass rate ON, Munk moment ON,
% LEV OFF, application point 'colocated', Rossby gate 'kinematic'.
% Switching rate + moment + edge drag OFF leaves the CORE shared with the planar
% model (APW coefficients + strip theory + rigid-body + added mass) -- see
% testing/shape3d/testShape3DFlatEquivalence.m.
% cfg.enableAddedMassRate   = false;   % Adot*v      (APW eqs. 6.1-6.2 cross terms)
% cfg.enableAddedMassMoment = false;   % v x (A v)   (APW eq. 6.3, the Munk moment)
% cfg.enableEdgeDrag        = false;   % tip crossflow drag
% cfg.enableLEV             = true;    % leading-edge-vortex vortex lift
% cfg.levApplicationPoint   = 'forward';
% cfg.lev = struct('rossbyDefinition','geometric');

% --- Build ----------------------------------------------------------------
bsp.seedShape     = polyshape([-hs, hs, hs, -hs], [-hc, -hc, hc, hc]);  % drawing x = span
bsp.seedDensity   = arealDensity;
bsp.seedThickness = thickness;
bsp.numStrips     = numStrips;
bsp.tSamples      = tSamples;
bsp.nutMass_t     = nutMass * ones(size(tSamples));
bsp.nutPos_t      = repmat(nutPos, 1, numel(tSamples));
if ~isempty(twistFn);     bsp.twist     = twistFn;     end
if ~isempty(curvatureFn); bsp.curvature = curvatureFn; end

seedParams = buildSeedParams(bsp, cfg);

fprintf('\nmodel: %s | strips %d | twist %g deg | curvature %g deg\n', ...
        seedParams.model, numStrips, tipTwistDeg, tipCurveDeg);
fprintf('physics: edgeDrag %d, addedMassRate %d, addedMassMoment %d, LEV %d (%s, %s)\n', ...
        seedParams.enableEdgeDrag, seedParams.enableAddedMassRate, ...
        seedParams.enableAddedMassMoment, seedParams.enableLEV, ...
        seedParams.levApplicationPoint, seedParams.lev.rossbyDefinition);
fprintf('mass: total %.3g kg, CoM [%+.4f %+.4f %+.4f] m\n', ...
        seedParams.massParams.M_total, seedParams.massParams.com_t(:,1));

%% 2. Visualize the seed shape (sanity check before integrating)
% Worth a look whenever twist or curvature is on: this is where a sign error in a
% shape profile shows up as an obviously wrong plate, rather than as a puzzling
% trajectory ten minutes later.
figure('Name', 'Seed shape sanity check (3D)');
visualizeSeedShape(seedParams, eye(3), [0;0;0], struct());
axis equal; grid on;
xlim([-0.05 0.05]); ylim([-0.05 0.05]); zlim([-0.05 0.05]);
xlabel('World X (m)'); ylabel('World Z (m)'); zlabel('World Y (m) -- up');
title(sprintf('Seed geometry, identity pose  (twist %g deg, curvature %g deg)', ...
              tipTwistDeg, tipCurveDeg));
view(35, 20);

%% 3. Run the dynamics (ode45)
r0     = [0; 0; 0];
q0     = [1; 0; 0; 0];              % level release. Tilted: axisAngleToQuat([0;0;1], pi/6)
v0     = [0; 0; 0];
omega0 = [0; 10; 0];                % initial spin about the plate normal (rad/s)

x0    = [r0; q0; v0; omega0];
tspan = 0:0.003:10;                 % integration horizon (s)

odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
rhs     = seedRHS(seedParams);      % -> seed6DOFODE3D for a shape3d seed
[tOut, xOut] = ode45(@(t, x) rhs(t, x, seedParams), tspan, x0, odeOpts);

% xOut columns: [ r(1:3), q(4:7), v(8:10), omega(11:13) ]. ode45 discards the
% optional `intermediates` output of the RHS -- if per-strip quantities (alpha,
% CT, CD, the LEV gate, tau_addedMass, ...) are needed, re-call the RHS on the
% rows of xOut in post-processing.
fprintf('\nintegrated %d steps, final height %+.3f m\n', numel(tOut), xOut(end,2));

%% 4. Visualize the trajectory, and print the metrics + label
visualizeSeedTrajectory(tOut, xOut(:,1:3).', xOut(:,4:7).');

% The same metrics and classifier the sweeps use, on this one run. Printed here
% so a case can be checked BY EYE against the label it would be given -- which is
% the point of having a single-drop harness at all.
mopts = struct('windowStartFrac', 0.5, 'convergeTol', 0.20);
m     = computeTrajectoryMetrics(tOut, xOut, mopts);
[lab, labInfo] = classifyFlightMode(m, defaultModeThresholds());
fprintf(['\nmode = %s\n  descent %.2f m/s | vSpin %.2f | spanSpin %.2f rad/s\n' ...
         '  cone %.0f deg | glide %.2f | converged %d\n'], ...
        lab, m.descentSpeed, m.verticalSpinMag, m.spanwiseSpin, ...
        m.coneAngleDeg, m.glideRatio, m.converged);
if isfield(labInfo, 'reason'); fprintf('  reason: %s\n', labInfo.reason); end

%% 5. Animate and save
% Set outFolder OUTSIDE the repo: generated videos stay out of git.
outFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Compiled_Outputs";
if ~exist(outFolder, 'dir'); mkdir(outFolder); end
runName = 'quick_Spin_3D';

traj.t             = tOut;
traj.quat          = xOut(:, 4:7).';
traj.comPos_world  = xOut(:, 1:3).';
traj.comVel_world  = xOut(:, 8:10).';
traj.comOmega_body = xOut(:, 11:13).';
% comPos_body omitted -> animateSeed falls back to massParams.com_t.

opts = struct();
opts.outputFolder = char(outFolder);
opts.videoFile    = char(fullfile(outFolder, [runName '.mp4']));
opts.fps          = 10;
opts.axisMode     = 'follow';
opts.zoomFactor   = 3;
opts.frameStep    = 5;
opts.showVels     = false;
opts.velOpts.showTotal     = false;
opts.velOpts.showProjected = false;
opts.velOpts.scale         = 0.04;
opts.figSize       = [600 400];
opts.trajLineWidth = 1;
opts.shapeOpts.showCom = false;
opts.shapeOpts.showNut = false;

animateSeed(seedParams, traj, opts);

seedRunData = struct('seedParams',seedParams, 'x0',x0, 'tOut',tOut, 'xOut',xOut, ...
                     'cfg',cfg, 'mode',lab, 'metrics',m, 'currDate',datetime('now'));
save(fullfile(outFolder, [runName '.mat']), 'seedRunData');
fprintf('\nwrote %s.mp4 / .mat to %s\n', runName, outFolder);
