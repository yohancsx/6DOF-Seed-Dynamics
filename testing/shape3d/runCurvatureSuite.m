%% Curvature test suite (3D-shape model) -- what does spanwise curvature do?
% Companion to runTwistSuite, exercising the CURVATURE path of the non-planar
% model (setupSeedShape3D + seed6DOFODE3D, reached through the shared
% buildSeedParams / runSingleMode via the model tag). The nut is at the CENTRE of
% the plate in every case, so anything that happens comes purely from the SHAPE.
%
% Curvature is dihedral: a rotation of each strip about the CHORD axis, so unlike
% twist the strip centres LEAVE the body x-z plane. The flat span coordinate is
% treated as ARC LENGTH and the tangent integrated, so the wing keeps its true
% length as it bends. Expect curvature to act mainly on the DESCENT / cone (a bowl
% presents more area downward and moves the CoM out of plane) rather than on spin,
% which is twist's mechanism.
%
% Two families, each swept small -> large tip dihedral:
%   A. SYMMETRIC BOWL: phi(s) = k*s -- both half-wings tilt the same way in y,
%                      giving a symmetric bowl (y(s) is even).
%   B. SINGLE SIDE   : only the +s half-wing bends, phi(s) = (s>0)*k*s; the -s
%                      half stays flat.
%
% ---------------------------------------------------------------------------
% THE FLAT REFERENCE IS NOW A SANITY CHECK. This suite used to carry a caveat: a
% zero-valued curvature profile took the curved branch, which switched the old
% span force off, so the 0 deg rows ran different physics from a genuinely flat
% seed. Phase 4 retired the span force from the 3D model, so that confound is
% gone -- the 0 deg rows and the flat reference now run IDENTICAL physics and
% must agree (to integrator precision; the curved-mass path reproduces the flat
% mass properties to ~1e-14, see testCurvedMassProperties check B). Section 4
% prints the mismatch; anything beyond noise is a bug.
% ---------------------------------------------------------------------------
%
% Section-organised (no local functions), so it saves as an .mlx cleanly.
% Body axes: x = chord, y = normal, z = span.

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
% (No cfg.aero: the shape3d model has no span force or span torque left for
%  C_span / C_span_torque to scale. Its physics defaults -- edge drag and the
%  added-mass rate ON -- come from setupSeedShape3D.)

% --- Simulation / analysis ------------------------------------------------
cfg.tspan = [0 12];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();

% --- Curvature sweep + release --------------------------------------------
tipDeg    = [0 5 10 15 20 25 30 35];   % tip dihedral magnitudes (deg)
q0        = [1; 0; 0; 0];              % released level
omega0    = [0; 0; 0];                 % from rest
nutCentre = [0; 0; 0];                 % CoM at the plate centre (no mass offset)
showTraj  = true;                      % show a couple of representative trajectories

% --- Base-seed-params -----------------------------------------------------
xh = cfg.spanLength/2;   yh = cfg.chordLength/2;   hs = cfg.spanLength/2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

%% 1. Flat reference -- no curvature field at all (the 0 deg rows must match it)
fprintf('\n==== Flat reference (no curvature field) ====\n');
rFlat    = runSingleMode('flat_reference', nutCentre, q0, omega0, cfg, baseBsp);
comYflat = rFlat.seedParams.massParams.com_t(2,1);
fprintf('  CoM_y = %+.3e m\n', comYflat);

%% 2. Family A -- symmetric bowl, phi(s) = k*s
fprintf('\n==== Family A: SYMMETRIC BOWL (phi = k*s) ====\n');
vSpinA = zeros(size(tipDeg));  spanA = zeros(size(tipDeg));  coneA = zeros(size(tipDeg));
descA  = zeros(size(tipDeg));  tiltA = zeros(size(tipDeg));  comYA = zeros(size(tipDeg));
modeA  = strings(size(tipDeg));
resA   = struct('t',{},'x',{},'sp',{});
for j = 1:numel(tipDeg)
    td = tipDeg(j);   k = deg2rad(td)/hs;
    bb = baseBsp;   bb.curvature = @(s) k .* s;          % phi(s) = k*s
    r  = runSingleMode(sprintf('A_bowl_%02ddeg', td), nutCentre, q0, omega0, cfg, bb);
    vSpinA(j)=r.metrics.verticalSpinMag; spanA(j)=r.metrics.spanwiseSpin;
    coneA(j)=r.metrics.coneAngleDeg;     descA(j)=r.metrics.descentSpeed;
    tiltA(j)=r.metrics.tiltStd;          modeA(j)=string(r.mode);
    comYA(j)=r.seedParams.massParams.com_t(2,1);
    resA(j).t=r.t; resA(j).x=r.x; resA(j).sp=r.seedParams;
end

%% 3. Family B -- single side only (+s half bent)
fprintf('\n==== Family B: SINGLE SIDE (+s half-wing bent) ====\n');
vSpinB = zeros(size(tipDeg));  spanB = zeros(size(tipDeg));  coneB = zeros(size(tipDeg));
descB  = zeros(size(tipDeg));  tiltB = zeros(size(tipDeg));  comYB = zeros(size(tipDeg));
modeB  = strings(size(tipDeg));
resB   = struct('t',{},'x',{},'sp',{});
for j = 1:numel(tipDeg)
    td = tipDeg(j);   k = deg2rad(td)/hs;
    bb = baseBsp;   bb.curvature = @(s) (s > 0) .* k .* s;   % only +s half bent
    r  = runSingleMode(sprintf('B_oneSide_%02ddeg', td), nutCentre, q0, omega0, cfg, bb);
    vSpinB(j)=r.metrics.verticalSpinMag; spanB(j)=r.metrics.spanwiseSpin;
    coneB(j)=r.metrics.coneAngleDeg;     descB(j)=r.metrics.descentSpeed;
    tiltB(j)=r.metrics.tiltStd;          modeB(j)=string(r.mode);
    comYB(j)=r.seedParams.massParams.com_t(2,1);
    resB(j).t=r.t; resB(j).x=r.x; resB(j).sp=r.seedParams;
end

%% 4. Summary tables
fprintf('\n--- Family A: symmetric bowl ---\n');
fprintf('%-5s %-9s %-9s %-8s %-7s %-12s\n', ...
        'tip','CoM_y (m)','vSpin','descent','cone','mode');
for j = 1:numel(tipDeg)
    fprintf('%-5d %-+9.2e %-9.2f %-8.2f %-7.1f %-12s\n', ...
            tipDeg(j), comYA(j), vSpinA(j), descA(j), coneA(j), modeA(j));
end

fprintf('\n--- Family B: single side ---\n');
fprintf('%-5s %-9s %-9s %-8s %-7s %-12s\n', ...
        'tip','CoM_y (m)','vSpin','descent','cone','mode');
for j = 1:numel(tipDeg)
    fprintf('%-5d %-+9.2e %-9.2f %-8.2f %-7.1f %-12s\n', ...
            tipDeg(j), comYB(j), vSpinB(j), descB(j), coneB(j), modeB(j));
end

fprintf('\n--- Flat reference ---\n');
fprintf('%-5s %-+9.2e %-9.2f %-8.2f %-7.1f %-12s\n', ...
        'flat', comYflat, rFlat.metrics.verticalSpinMag, ...
        rFlat.metrics.descentSpeed, rFlat.metrics.coneAngleDeg, string(rFlat.mode));

% Sanity check: the 0 deg rows run identical physics to the flat reference.
i0 = find(tipDeg == 0, 1);
if ~isempty(i0)
    dA = abs(descA(i0) - rFlat.metrics.descentSpeed);
    dB = abs(descB(i0) - rFlat.metrics.descentSpeed);
    fprintf(['\nSanity: |descent(0 deg) - descent(flat)| = %.2e (bowl), %.2e (single side)' ...
             '  -- should be integrator noise.\n'], dA, dB);
end

%% 5. Plots
figure('Name','Curvature sweep','Color','w');

subplot(2,2,1);
plot(tipDeg, descA, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, descB, '-s', 'LineWidth',1.5);
yline(rFlat.metrics.descentSpeed, '--', 'flat ref', 'Color',[0.45 0.45 0.45]);
grid on; xlabel('tip dihedral (deg)'); ylabel('descent speed (m/s)');
legend('symmetric bowl','single side','Location','best'); title('Descent speed vs curvature');

subplot(2,2,2);
plot(tipDeg, vSpinA, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, vSpinB, '-s', 'LineWidth',1.5);
yline(rFlat.metrics.verticalSpinMag, '--', 'flat ref', 'Color',[0.45 0.45 0.45]);
grid on; xlabel('tip dihedral (deg)'); ylabel('vertical spin |\omega_Y^{world}| (rad/s)');
legend('symmetric bowl','single side','Location','best'); title('Vertical-axis spin vs curvature');

subplot(2,2,3);
plot(tipDeg, coneA, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, coneB, '-s', 'LineWidth',1.5);
grid on; xlabel('tip dihedral (deg)'); ylabel('cone angle (deg)');
legend('symmetric bowl','single side','Location','best'); title('Cone angle vs curvature');

subplot(2,2,4);
plot(tipDeg, comYA*1e3, '-o', 'LineWidth',1.5); hold on;
plot(tipDeg, comYB*1e3, '-s', 'LineWidth',1.5);
grid on; xlabel('tip dihedral (deg)'); ylabel('CoM offset from plane, y (mm)');
legend('symmetric bowl','single side','Location','best'); title('Out-of-plane CoM shift');

%% 6. Representative trajectories (largest curvature of each family)
if showTraj
    jA = numel(tipDeg);
    visualizeSeedTrajectory(resA(jA).t, resA(jA).x(:,1:3).', resA(jA).x(:,4:7).', ...
        struct('fig1Name',sprintf('A bowl %ddeg - traj',tipDeg(jA)), ...
               'fig2Name',sprintf('A bowl %ddeg - angles',tipDeg(jA))));
    visualizeSeedTrajectory(resB(jA).t, resB(jA).x(:,1:3).', resB(jA).x(:,4:7).', ...
        struct('fig1Name',sprintf('B single-side %ddeg - traj',tipDeg(jA)), ...
               'fig2Name',sprintf('B single-side %ddeg - angles',tipDeg(jA))));
end

fprintf(['\nDone. (Expect the CoM to move out of plane with curvature -- that is the\n' ...
         'signature twist does NOT have -- and the descent/cone to respond more strongly\n' ...
         'than the spin. Mass is conserved throughout; see testCurvedMassProperties.)\n']);
