%% Span-force comparison -- 4 hand-tuned cases, span forces OFF vs ON
% Isolates the span-force question with four single runs you can eyeball:
%
%   Case 1  span OFF  tumbling      Case 3  span ON   tumbling
%   Case 2  span OFF  autorotation  Case 4  span ON   autorotation
%
% Each case runs the dynamics and calls visualizeSeedTrajectory, so you get the
% 3D trajectory AND the Euler-angle history for every run.
%
% The seed geometry / mass / environment are built ONCE in Section 0 and reused;
% each case edits only the NUT POSITION (and the span-force switches). Tune the
% nut positions and the switches per case, then re-run that section.
%
% Body axes: x = chord, y = plate-normal (vertical when level), z = span.
% Nut position is given in HALF-DIMENSION fractions and converted to metres in
% each case, so the numbers stay comparable to the mode-suite sweeps:
%   fChord -> multiples of half-chord (hc), fSpan -> multiples of half-span (hs).
%   Values > 1 place the nut OFF the planform (allowed, and often needed).
%
% Free of local functions, so it can be saved as a live script (.mlx).

%% 0. Shared setup  -- geometry, environment, and the baseline seed
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ are assumed already on the MATLAB path.

% --- Baseline seed geometry / material ------------------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   % kg

% --- Environment ----------------------------------------------------------
cfg.rhoFluid = 1.225;   cfg.g = 9.81;

% --- Simulation -----------------------------------------------------------
cfg.tspan = [0 10];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;

% --- Span-force tuning knobs (shared defaults; override per case) ---------
% C_span       scales the span FORCE; C_span_torque scales the span TORQUE only.
% enableSpanCOPMigration = false puts the span torque at the fixed geometric
% centre instead of the migrating span CoP (removes the oscillating arm).
cfg.aero = struct('C_span', 1.0, 'C_span_torque', 1.0);   % partial override; rest = defaults
cfg.enableSpanGeomVelocity = true;
cfg.enableSpanCOPMigration = true;

% --- Half-dimension references (nut positions are given in these units) ---
hc = cfg.chordLength / 2;    % half-chord (body x)
hs = cfg.spanLength  / 2;    % half-span  (body z)

% --- Baseline base-seed-params (nut position overridden per case) ---------
xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;
baseBsp.tSamples      = cfg.tSamples;
baseBsp.nutMass_t     = cfg.nutMass * ones(size(cfg.tSamples));
baseBsp.nutPos_t      = repmat([0;0;0], 1, numel(cfg.tSamples));

odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);

%% 1. Span forces OFF  --  TUMBLING / spiralling
% ---- EDIT: nut position (fractions of half-chord / half-span) ------------
fChord = 0.0;    % chordwise offset
fSpan  = 0.4;    % spanwise offset  (modest -> spiral tumbling regime)
% -------------------------------------------------------------------------
cfgCase = cfg;   cfgCase.enableSpanForce = false;
bsp = baseBsp;   bsp.nutPos_t = repmat([0; 0; cfg.spanLength*0.05], 1, numel(cfg.tSamples));
sp  = buildSeedParams(bsp, cfgCase);

q0 = axisAngleToQuat([0; 0; 1], pi/6);      % pi/6 tilt about body z, to start the tumble
x0 = [zeros(3,1); q0; zeros(3,1); [0;0;0]];

[t1, x1] = ode45(@(t,x) seed6DOFODE(t, x, sp), cfg.tspan, x0, odeOpts);
visualizeSeedTrajectory(t1, x1(:,1:3).', x1(:,4:7).', ...
    struct('fig1Name','C1 spanOFF tumbling - traj', 'fig2Name','C1 spanOFF tumbling - angles'));
fprintf('Case 1 (span OFF, tumbling): %d steps, final height %.3f m\n', numel(t1), x1(end,2));

%% 2. Span forces OFF  --  AUTOROTATION
% ---- EDIT: nut position (fractions of half-chord / half-span) ------------
fChord = 0.87;   % chordwise offset (~puts CoM near the chordwise CoP)
fSpan  = 1.5;    % spanwise offset  (tip-heavy -> autorotation-prone)
% -------------------------------------------------------------------------
cfgCase = cfg;   cfgCase.enableSpanForce = false;
bsp = baseBsp;   bsp.nutPos_t = repmat([cfg.chordLength; 0; cfg.spanLength], 1, numel(cfg.tSamples));
sp  = buildSeedParams(bsp, cfgCase);

q0 = axisAngleToQuat([0; 0; 1], 0);
x0 = [zeros(3,1); q0; zeros(3,1); [0; 2; 0]];   % small spin nudge about the normal

[t2, x2] = ode45(@(t,x) seed6DOFODE(t, x, sp), cfg.tspan, x0, odeOpts);
visualizeSeedTrajectory(t2, x2(:,1:3).', x2(:,4:7).', ...
    struct('fig1Name','C2 spanOFF autorotation - traj', 'fig2Name','C2 spanOFF autorotation - angles'));
fprintf('Case 2 (span OFF, autorotation): %d steps, final height %.3f m\n', numel(t2), x2(end,2));

%% 3. Span forces ON  --  TUMBLING / spiralling
% ---- EDIT: nut position (fractions of half-chord / half-span) ------------
% -------------------------------------------------------------------------
cfgCase = cfg;   
cfgCase.enableSpanForce = true;
cfgCase.enableSpanCOPMigration = true;   % <-- set false to kill the migrating arm
cfgCase.aero.C_span            = 0.2;    % <-- span FORCE scale 0.2
cfgCase.aero.C_span_torque     = 0.7;    % <-- span TORQUE scale (independent) 0.7

bsp = baseBsp;   bsp.nutPos_t = repmat([0; 0; cfg.spanLength*0.05], 1, numel(cfg.tSamples));
sp  = buildSeedParams(bsp, cfgCase);

q0 = axisAngleToQuat([0; 0; 1], pi/6);
x0 = [zeros(3,1); q0; zeros(3,1); [0;0;0]];

[t3, x3] = ode45(@(t,x) seed6DOFODE(t, x, sp), cfg.tspan, x0, odeOpts);
visualizeSeedTrajectory(t3, x3(:,1:3).', x3(:,4:7).', ...
    struct('fig1Name','C3 spanON tumbling - traj', 'fig2Name','C3 spanON tumbling - angles'));
fprintf('Case 3 (span ON, tumbling): %d steps, final height %.3f m\n', numel(t3), x3(end,2));

%% 4. Span forces ON  --  AUTOROTATION
% This is the case that "almost gets there" then destabilises. Try, in order:
%   (a) cfgCase.enableSpanCOPMigration = false   -> fixed geometric-centre arm
%   (b) cfgCase.aero.C_span_torque = 0.5 / 0.2   -> weaken the torque only
%   (c) cfgCase.aero.C_span        = 0.5         -> weaken the force too
% ---- EDIT: nut position (fractions of half-chord / half-span) ------------
% ---- EDIT: span-force switches for this case ----------------------------
cfgCase = cfg;
cfgCase.enableSpanForce        = true;
cfgCase.enableSpanCOPMigration = false;   % <-- set false to kill the migrating arm
cfgCase.aero.C_span            = 0.2;    % <-- span FORCE scale 0.2
cfgCase.aero.C_span_torque     = 0.7;    % <-- span TORQUE scale (independent) 0.7
% -------------------------------------------------------------------------
bsp = baseBsp;   
bsp.nutPos_t = repmat([cfg.chordLength; 0; cfg.spanLength], 1, numel(cfg.tSamples));
sp  = buildSeedParams(bsp, cfgCase);

q0 = axisAngleToQuat([0; 0; 1], 0);
x0 = [zeros(3,1); q0; zeros(3,1); [0; 0; 0]];

[t4, x4] = ode45(@(t,x) seed6DOFODE(t, x, sp), cfg.tspan, x0, odeOpts);
visualizeSeedTrajectory(t4, x4(:,1:3).', x4(:,4:7).', ...
    struct('fig1Name','C4 spanON autorotation - traj', 'fig2Name','C4 spanON autorotation - angles'));
fprintf('Case 4 (span ON, autorotation): %d steps, final height %.3f m\n', numel(t4), x4(end,2));

%% 5. TORQUE-BUDGET DIAGNOSTIC  --  do the strips produce roll on their own?
% Self-contained: builds ONE seed, integrates it, then re-calls seed6DOFODE on
% every saved state to recover the (ode45-discarded) intermediate torques, and
% splits each body-axis torque into its STRIP-sum vs WHOLE-SEED-SPAN parts.
%
% The question this answers: during tumbling/spiral, is the roll torque (body x)
% coming mainly from the STRIPS (their spanwise-distributed normal loading) or
% from the SPAN-force torque? If the strips already supply the roll, the span
% torque is redundant and can be dropped. Compare a tumbling and an
% autorotation config by editing the block below.
%
% Body-axis torque components:  tau(1)=roll (x)  tau(2)=yaw (y)  tau(3)=pitch (z)

% ---- EDIT: which configuration to dissect --------------------------------
cfgD = cfg;
cfgD.enableSpanForce        = true;
cfgD.enableSpanCOPMigration = false;
cfgD.aero.C_span            = 0.2;
cfgD.aero.C_span_torque     = 0.7;


nutPosD = [0; 0; cfg.spanLength*0.05];     % tumbling-prone nut position
q0D     = axisAngleToQuat([0; 0; 1], pi/6);
omega0D = [0; 0; 0];

nutPosD = [cfg.chordLength; 0; cfg.spanLength];     % autorotation
q0D     = axisAngleToQuat([0; 0; 1], 0);
omega0D = [0; 0; 0];
% -------------------------------------------------------------------------

bspD = baseBsp;   bspD.nutPos_t = repmat(nutPosD, 1, numel(cfg.tSamples));
spD  = buildSeedParams(bspD, cfgD);
x0D  = [zeros(3,1); q0D; zeros(3,1); omega0D];
[tD, xD] = ode45(@(t,x) seed6DOFODE(t, x, spD), cfg.tspan, x0D, odeOpts);

% --- Recover per-state torque contributions (ode45 dropped the 2nd output) -
nD       = numel(tD);
tauStrip = zeros(3, nD);   % strip-sum torque (body)
tauSpan  = zeros(3, nD);   % whole-seed span-force torque (body)
tauTot   = zeros(3, nD);   % applied total torque (body)
rSpanCoP = zeros(3, nD);   % span-CoP application point rel. to CoM (body frame)
FspanFull= zeros(3, nD);   % full span force [0; Fy; Fz] (body frame)
for k = 1:nD
    [~, im] = seed6DOFODE(tD(k), xD(k, :).', spD);
    tauStrip(:, k)  = sum(im.tau_strip_body, 2);
    tauSpan(:, k)   = im.tau_span;
    tauTot(:, k)    = im.tau_body;
    rSpanCoP(:, k)  = im.r_spanCoP_body;
    FspanFull(:, k) = im.F_span_full;
end

% --- Plot: strip vs span vs total, one panel per body axis ----------------
axisNames = {'Roll (body x)', 'Yaw (body y)', 'Pitch (body z)'};
figure('Name', 'Torque budget: strips vs span force');
for a = 1:3
    subplot(3, 1, a);
    plot(tD, tauStrip(a, :), 'LineWidth', 1.5); hold on;
    plot(tD, tauSpan(a, :),  'LineWidth', 1.5);
    plot(tD, tauTot(a, :),   'k--', 'LineWidth', 1.0);
    grid on; ylabel('N\cdotm');
    title(axisNames{a}, 'Interpreter', 'tex');
    if a == 1
        legend('strips (sum)', 'span force', 'applied total', 'Location', 'best');
    end
end
xlabel('Time (s)');

% --- Steady-state roll budget (latter half): a single number to read ------
win   = tD >= 0.5 * tD(end);
rmsSt = rms(tauStrip(1, win));
rmsSp = rms(tauSpan(1, win));
fprintf(['\nRoll-torque budget (latter half, RMS):\n', ...
         '  strips = %.3e N*m\n  span   = %.3e N*m\n  ratio span/strips = %.2f\n'], ...
         rmsSt, rmsSp, rmsSp / max(rmsSt, eps));
if rmsSt > rmsSp
    fprintf('  -> STRIPS dominate roll: the span torque is likely redundant.\n');
else
    fprintf('  -> SPAN torque dominates roll: strips under-produce it.\n');
end

% =========================================================================
% ANIMATION: seed body + CoM + spanwise CoP + F_span_full normal (y) direction
% =========================================================================
% Debug view for the span-CoP hypothesis. Each frame draws the seed at its pose
% and overlays:
%   red dot     -- centre of mass (CoM)
%   blue dot    -- spanwise centre of pressure (where the span torque is applied)
%   magenta arrow at the CoP -- direction of F_span_full's NORMAL (body-y)
%                  component (the component that produces the roll torque)
% Watch whether the CoP sits where it should relative to the seed (span centre?
% migrating to a tip? oscillating?), and which way the normal force points.
%
% Same world->plot remap ([X; Z; Y]) as the visualizers, so overlays align with
% the seed body drawn by visualizeSeedShape.
toPlot  = @(w) [w(1,:); w(3,:); w(2,:)];
comBody = spD.massParams.com_t(:, 1);                 % CoM in body datum (nut fixed here)
longestDim = max(cfg.spanLength, cfg.chordLength);
halfWin    = 2.5 * longestDim;                         % follow-cube half-window

% Scale the normal-force arrow so its peak spans ~half the seed span (direction
% and RELATIVE magnitude are what matter, not the absolute Newtons).
FyMax   = max(abs(FspanFull(2, :)));
arrowScale = (0.5 * cfg.spanLength) / max(FyMax, eps);

frameStride = max(1, round(nD / 150));                 % ~150 frames
shapeOpts = struct('showCom', false, 'showNut', false, 'showGeoCenter', false, ...
                   'showPlate', true, 'showPlateEdges', true);

figAnim = figure('Name', 'Span CoP / CoM / F_span,y animation', 'Color', 'w');
ax = axes('Parent', figAnim); hold(ax, 'on');
for k = 1:frameStride:nD
    Rk       = quatToRotm(xD(k, 4:7).');
    comW     = xD(k, 1:3).';                            % CoM in world
    posDatum = comW - Rk * comBody;                     % body-datum origin in world

    cla(ax); hold(ax, 'on');
    visualizeSeedShape(spD, Rk, posDatum, shapeOpts);   % seed body (draws into current axes)

    % CoM (red) and span CoP (blue), both mapped to plot coords.
    comP = toPlot(comW);
    copW = comW + Rk * rSpanCoP(:, k);                  % span-CoP world position
    copP = toPlot(copW);
    hCom = plot3(comP(1), comP(2), comP(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    hCop = plot3(copP(1), copP(2), copP(3), 'bd', 'MarkerFaceColor', 'b', 'MarkerSize', 9);

    % F_span_full normal (body-y) component as an arrow at the CoP (magenta).
    Fdir_world = Rk * [0; FspanFull(2, k) * arrowScale; 0];
    Fdir_plot  = toPlot(Fdir_world);
    hF = quiver3(copP(1), copP(2), copP(3), Fdir_plot(1), Fdir_plot(2), Fdir_plot(3), ...
                 0, 'm', 'LineWidth', 2, 'MaxHeadSize', 1);

    % Follow the seed with an equal-aspect cube so proportions stay true.
    xlim(ax, comP(1) + [-halfWin halfWin]);
    ylim(ax, comP(2) + [-halfWin halfWin]);
    zlim(ax, comP(3) + [-halfWin halfWin]);
    daspect(ax, [1 1 1]); grid(ax, 'on'); view(ax, 35, 20);
    xlabel(ax, 'World X (m)'); ylabel(ax, 'World Z (m)'); zlabel(ax, 'World Y (m) -- up');
    title(ax, sprintf('t = %.2f s   |   span-CoP z-offset from CoM = %+.4f m', ...
                      tD(k), rSpanCoP(3, k)), 'Interpreter', 'none');
    drawnow;
end
legend([hCom, hCop, hF], {'CoM', 'span CoP', 'F_{span} normal (y)'}, 'Location', 'best');
