%% Twist test suite (3D-shape model) -- does spanwise twist spin a centred seed?
% Exercises the twist path of the non-planar model (setupSeedShape3D +
% seed6DOFODE3D, reached through the shared buildSeedParams / runSingleMode via
% the model tag). The nut is at the CENTRE of the plate in every case, so any
% vertical-axis spin comes purely from the SHAPE, not a mass offset -- the samara
% autorotation mechanism.
%
% Two families, each swept small -> large tip angle:
%   A. BOTH ENDS  : anti-symmetric linear twist theta(z) = k*z (one half-wing
%                   pitched up, the other down -- a propeller twist).
%   B. SINGLE SIDE: only the +z half-wing twisted, theta(z) = (z>0)*k*z; the -z
%                   half stays flat.
%
% For each case it prints a metric/mode summary and collects the vertical spin,
% then plots spin vs tip twist for both families. Body axes: x=chord, y=normal,
% z=span. Section-organised (no local functions), so it saves as an .mlx cleanly.

%% 0. Configuration  -- EDIT HERE
root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% --- Base seed (the working seed) -----------------------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   cfg.rhoFluid = 1.225;      cfg.g = 9.81;
cfg.shapeModel  = 'shape3d';                          % <-- selects the 3D model
cfg.aero        = struct('C_span', 0.2, 'C_span_torque', 0.7);

% --- Simulation / analysis ------------------------------------------------
cfg.tspan = [0 12];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();

% --- Twist sweep + release ------------------------------------------------
tipDeg    = [0 5 10 15 20 25 30];   % tip twist magnitudes (deg)
q0        = [1; 0; 0; 0];           % released level
omega0    = [0; 0; 0];              % from rest
nutCentre = [0; 0; 0];              % CoM at the plate centre (no mass offset)
showTraj  = true;                   % show a couple of representative trajectories

% --- Base-seed-params -----------------------------------------------------
xh = cfg.spanLength/2;   yh = cfg.chordLength/2;   hs = cfg.spanLength/2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

%% 1. Family A -- both ends (anti-symmetric linear twist)
fprintf('\n==== Family A: BOTH ENDS (anti-symmetric linear twist) ====\n');
vSpinA = zeros(size(tipDeg));  spanA = zeros(size(tipDeg));  coneA = zeros(size(tipDeg));
tiltA  = zeros(size(tipDeg));  modeA = strings(size(tipDeg));  resA = struct('t',{},'x',{},'sp',{});
for j = 1:numel(tipDeg)
    td = tipDeg(j);   k = deg2rad(td)/hs;
    bb = baseBsp;   bb.twist = @(z) k*z;                 % theta(z) = k*z
    r  = runSingleMode(sprintf('A_bothEnds_%02ddeg', td), nutCentre, q0, omega0, cfg, bb);
    vSpinA(j)=r.metrics.verticalSpinMag; spanA(j)=r.metrics.spanwiseSpin;
    coneA(j)=r.metrics.coneAngleDeg; tiltA(j)=r.metrics.tiltStd; modeA(j)=string(r.mode);
    resA(j).t=r.t; resA(j).x=r.x; resA(j).sp=r.seedParams;
end

%% 2. Family B -- single side only (+z half twisted)
fprintf('\n==== Family B: SINGLE SIDE (+z half-wing twisted) ====\n');
vSpinB = zeros(size(tipDeg));  spanB = zeros(size(tipDeg));  coneB = zeros(size(tipDeg));
tiltB  = zeros(size(tipDeg));  modeB = strings(size(tipDeg));  resB = struct('t',{},'x',{},'sp',{});
for j = 1:numel(tipDeg)
    td = tipDeg(j);   k = deg2rad(td)/hs;
    bb = baseBsp;   bb.twist = @(z) (z > 0) .* k .* z;   % only +z half twisted
    r  = runSingleMode(sprintf('B_oneSide_%02ddeg', td), nutCentre, q0, omega0, cfg, bb);
    vSpinB(j)=r.metrics.verticalSpinMag; spanB(j)=r.metrics.spanwiseSpin;
    coneB(j)=r.metrics.coneAngleDeg; tiltB(j)=r.metrics.tiltStd; modeB(j)=string(r.mode);
    resB(j).t=r.t; resB(j).x=r.x; resB(j).sp=r.seedParams;
end

%% 3. Summary tables + plots
fprintf('\n%-6s | %-22s | %-22s\n','tip','Family A (both ends)','Family B (single side)');
fprintf('%-6s | %-10s %-11s | %-10s %-11s\n','deg','vSpin','mode','vSpin','mode');
for j = 1:numel(tipDeg)
    fprintf('%-6d | %-10.1f %-11s | %-10.1f %-11s\n', tipDeg(j), vSpinA(j), modeA(j), vSpinB(j), modeB(j));
end

figure('Name','Twist-induced spin','Color','w');
subplot(1,2,1);
plot(tipDeg, vSpinA, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, vSpinB, '-s', 'LineWidth',1.5);
grid on; xlabel('tip twist (deg)'); ylabel('vertical spin  |\omega_Y^{world}|  (rad/s)');
legend('both ends','single side','Location','northwest'); title('Vertical-axis spin vs twist');
subplot(1,2,2);
plot(tipDeg, spanA, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, spanB, '-s', 'LineWidth',1.5);
grid on; xlabel('tip twist (deg)'); ylabel('spanwise spin  |\omega_z^{body}|  (rad/s)');
legend('both ends','single side','Location','best'); title('Spanwise (tumbling) spin vs twist');

%% 4. Representative trajectories (largest twist of each family)
if showTraj
    jA = numel(tipDeg);
    visualizeSeedTrajectory(resA(jA).t, resA(jA).x(:,1:3).', resA(jA).x(:,4:7).', ...
        struct('fig1Name',sprintf('A both-ends %ddeg - traj',tipDeg(jA)), ...
               'fig2Name',sprintf('A both-ends %ddeg - angles',tipDeg(jA))));
    visualizeSeedTrajectory(resB(jA).t, resB(jA).x(:,1:3).', resB(jA).x(:,4:7).', ...
        struct('fig1Name',sprintf('B single-side %ddeg - traj',tipDeg(jA)), ...
               'fig2Name',sprintf('B single-side %ddeg - angles',tipDeg(jA))));
end

fprintf('\nDone. (tip=0 must show ~0 vertical spin; spin should grow with twist.)\n');
