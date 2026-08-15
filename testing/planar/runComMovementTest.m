%% Time-varying CoM (nut-position) test  -- three scripted CoM motions
% Exercises the time-varying mass machinery (nutPos_t -> com_t, I_G_t, I_G_dot_t
% in setupSeedShapeAndMass, interpolated by getMassProperties, with the I_G_dot*w
% term active in rotationDynamics). Every earlier test used a single time sample,
% so this is the first real workout of that path.
%
% Scenarios (each is one drop; the NUT slides within the body as it falls):
%   1. Spanwise sweep : nut  0 -> +span -> 0 -> -span -> 0, dwelling at each.
%   2. Chordwise sweep: nut  0 -> +chord -> 0 -> -chord -> 0, dwelling at each.
%   3. Glide->spin->dive: nut starts chord-offset (glide), returns to centre,
%      slides spanwise (spin up), then chord-offset again (dive).
%
% Each scenario shows the 3D trajectory, the Euler angles, and the COMMANDED nut
% path vs time (so you can see the input that drove the motion).
%
% Edit the CoM limits (fractions of half-chord / half-span) and the dwell/move
% times in Section 0. Waypoints are interpolated with pchip (shape-preserving:
% flat dwells stay flat, moves ease in/out, no overshoot past the limits).
%
% CAVEAT: sliding an internal mass on a free-flying body is not momentum-rigorous
% here -- the model prescribes the CoM location within the body and tracks the
% CoM in inertial space, but does not add the reaction of the internal mass
% motion. This is a test of the code path and the mode transitions it drives, not
% a physically exact internal-actuation model.
%
% Uses local functions (at the end), so it stays self-contained; still converts
% to a live script.

%% 0. Configuration  -- EDIT HERE
helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ are assumed already on the MATLAB path.

% --- Baseline seed geometry / material ------------------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;
cfg.nutMass     = 75e-6;   % kg

% --- Environment + physics switches ---------------------------------------
cfg.rhoFluid = 1.225;   cfg.g = 9.81;
cfg.enableSpanForce            = true;
cfg.enableSpanTorqueAttenuation = true;
% (other enablers default as in buildSeedParams / setupSeedShapeAndMass)

% --- CoM travel limits (fractions of the half-dimension) ------------------
fChordMax  = 1.0;    % chordwise sweep amplitude (scenario 2), x half-chord
fSpanMax   = 0.05;    % spanwise  sweep amplitude (scenario 1), x half-span
% Scenario 3 waypoints:
fChordGlide = 0.8;   % initial chordwise offset -> glide
fSpanSpin   = 1.2;   % spanwise offset -> spin up
fChordDive  = 1.5;   % chordwise offset -> dive

% --- Timing ---------------------------------------------------------------
cfg.dwellTime = 1.0;    % pause held at each waypoint (s)
cfg.moveTime  = 0.8;    % transition time between waypoints (s)
cfg.dt        = 0.02;   % dense CoM-path / mass-sample timestep (s)

% --- Integrator -----------------------------------------------------------
cfg.odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% --- Combined mode animation (animateModeTrajectory) ----------------------
% Each scenario also renders a mode-coloured animation: 3D trajectory (left),
% CoM-within-the-seed-body (top-right), and a zoomed seed follow-cam showing the
% plate spinning/gliding (bottom-right). Set makeAnimation=false to skip the
% (slower) video writing. animFps sets the whole animation's frame rate;
% animPlaybackSpeed<1 slows it down (e.g. 0.25) to actually watch the fast spin.
cfg.makeAnimation     = true;
cfg.animFps           = 30;
cfg.animPlaybackSpeed = 0.25;
cfg.showSeedVels      = false;
here                  = fileparts(mfilename('fullpath'));
if isempty(here); here = pwd; end
cfg.animDir           = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\COM Movement Tests";   % where the .mp4s land

% --- Half-dimension references + baseline base-seed-params ----------------
hc = cfg.chordLength / 2;    hs = cfg.spanLength / 2;
xh = cfg.spanLength  / 2;    yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;
% (tSamples / nutMass_t / nutPos_t are filled per scenario by runComScenario)

%% 1. Spanwise CoM sweep  (0 -> +span -> 0 -> -span -> 0)
posList = [ [0;0;0], [0;0;+fSpanMax*hs], [0;0;0], [0;0;-fSpanMax*hs], [0;0;0] ];
dwell   = cfg.dwellTime * ones(1, size(posList,2));
runComScenario('Scenario 1: spanwise sweep', posList, dwell, cfg, baseBsp);

%% 2. Chordwise CoM sweep  (0 -> +chord -> 0 -> -chord -> 0)
posList = [ [0;0;0], [+fChordMax*hc;0;0], [0;0;0], [-fChordMax*hc;0;0], [0;0;0] ];
dwell   = cfg.dwellTime * ones(1, size(posList,2));
runComScenario('Scenario 2: chordwise sweep', posList, dwell, cfg, baseBsp);

%% 3. Glide -> centre -> spin -> dive
posList = [ [fChordGlide*hc; 0; 0], ...   % glide (chord offset)
            [0; 0; 0], ...                % back to centre
            [0; 0; fSpanSpin*hs], ...     % slide spanwise -> spin up
            [fChordDive*hc; 0; 0] ];      % chord offset -> dive
dwell   = cfg.dwellTime * ones(1, size(posList,2));
runComScenario('Scenario 3: glide-spin-dive', posList, dwell, cfg, baseBsp);


% =========================================================================
% LOCAL: build one CoM-movement scenario, integrate, and display it
% =========================================================================
function runComScenario(name, posList, dwellList, cfg, baseBsp)
% posList : 3xK nut-position waypoints, body [x;y;z] (m).
% dwellList: 1xK pause held at each waypoint (s). cfg.moveTime is the transition.

    % --- Smooth, dense nut-position path over the whole flight -------------
    [tD, nutPos, eventT] = buildComPath(posList, dwellList, cfg.moveTime, cfg.dt);

    % --- Build the time-varying seed --------------------------------------
    bsp = baseBsp;
    bsp.tSamples  = tD;
    bsp.nutPos_t  = nutPos;
    bsp.nutMass_t = cfg.nutMass * ones(size(tD));
    sp = buildSeedParams(bsp, cfg);

    % --- Integrate (from rest, level) -------------------------------------
    q0 = eul2quat([0 0 pi/6],'XYZ')';
    x0 = [zeros(3,1); q0; zeros(3,1); zeros(3,1)];
    [t, x] = ode45(@(tt,xx) seed6DOFODE(tt, xx, sp), [tD(1) tD(end)], x0, cfg.odeOpts);

    % --- Trajectory + Euler angles ----------------------------------------
    visualizeSeedTrajectory(t, x(:,1:3).', x(:,4:7).', ...
        struct('fig1Name', [name ' - traj'], 'fig2Name', [name ' - angles']));

    % --- Commanded nut path vs time (the input that drove the motion) -----
    figure('Name', [name ' - CoM command']);
    plot(tD, nutPos(1,:), 'LineWidth', 1.5); hold on;
    plot(tD, nutPos(2,:), 'LineWidth', 1.5);
    plot(tD, nutPos(3,:), 'LineWidth', 1.5);
    grid on; xlabel('Time (s)'); ylabel('Nut position (m)');
    legend('chord (body x)', 'normal (body y)', 'span (body z)', 'Location', 'best');
    title([name ' -- commanded nut position'], 'Interpreter', 'none');

    % --- Combined mode-coloured animation (3D traj + top view + Euler) -----
    % eventT (the arrival at each CoM waypoint) is marked, so the mode onsets
    % line up with the commanded moves.
    if isfield(cfg, 'makeAnimation') && cfg.makeAnimation
        vfile = fullfile(cfg.animDir, [matlab.lang.makeValidName(name) '.mp4']);
        animateModeTrajectory(t, x, struct( ...
            'videoFile',  vfile, ...
            'fps',        cfg.animFps, ...
            'playbackSpeed', cfg.animPlaybackSpeed, ...
            'showSeedVels', cfg.showSeedVels, ...
            'title',      name, ...
            'eventTimes', eventT, ...
            'seedParams', sp));
    end

    fprintf('%s: %.1f s flight, %d ode steps, final height %.3f m\n', ...
            name, tD(end), numel(t), x(end, 2));
end


% =========================================================================
% LOCAL: waypoint list (with dwells) -> dense, pchip-smoothed nut path
% =========================================================================
function [tDense, nutPos, eventT] = buildComPath(posList, dwellList, moveTime, dt)
% Builds a strictly-increasing waypoint schedule -- hold each posList column for
% its dwell, then move to the next over moveTime -- and interpolates it onto a
% dense grid with pchip (shape-preserving: flat dwells stay flat, no overshoot).
% eventT returns the arrival time of each waypoint (start of its dwell), for
% marking the CoM-move instants on downstream plots/animations.

    K   = size(posList, 2);
    wpT = 0;                    % waypoint times
    wpP = posList(:, 1);        % waypoint positions (3 x .)
    tcur = 0;
    eventT = zeros(1, K);       % arrival time of each waypoint (dwell start)
    for i = 1:K
        eventT(i) = tcur;                         % nut has arrived at waypoint i
        if dwellList(i) > 0                       % hold at waypoint i
            tcur = tcur + dwellList(i);
            wpT(end+1) = tcur;          %#ok<AGROW>
            wpP(:, end+1) = posList(:, i);         %#ok<AGROW>
        end
        if i < K                                  % move to waypoint i+1
            tcur = tcur + moveTime;
            wpT(end+1) = tcur;          %#ok<AGROW>
            wpP(:, end+1) = posList(:, i+1);       %#ok<AGROW>
        end
    end

    tDense = 0:dt:wpT(end);
    if tDense(end) < wpT(end); tDense(end+1) = wpT(end); end   % include the endpoint

    nutPos = zeros(3, numel(tDense));
    for a = 1:3
        nutPos(a, :) = interp1(wpT, wpP(a, :), tDense, 'pchip');
    end
end
