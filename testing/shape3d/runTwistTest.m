%% Twist animation test -- three twisted-seed drops, each rendered to an mp4
% Companion to runComMovementTest, but for the 3D-shape (twist) model. Each case
% drops a CENTRED-CoM seed with a different spanwise twist and renders the
% combined mode animation (3D trajectory + CoM-in-body + zoomed follow-cam that
% now shows the twisted strips) to a video, so you can regenerate the twisting
% animations in one run.
%
% Cases:
%   1. Slight twist        : anti-symmetric linear, small tip angle -> gentle spin.
%   2. Large twist         : anti-symmetric linear, large tip angle -> fast spin.
%   3. Small asymmetric     : only the +z half-wing twisted (single-side) -> spins
%      twist                 and drifts asymmetrically.
%
% The seed is built through the shared buildSeedParams (cfg.shapeModel='shape3d')
% and integrated via runSingleMode / seedRHS -> seed6DOFODE3D. Videos land in
% cfg.animDir (outside the repo). Section-organised with a local function, so it
% saves as an .mlx cleanly. Body axes: x=chord, y=normal, z=span.

%% 0. Configuration  -- EDIT HERE
root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

% --- Base seed (the working seed), 3D-shape model -------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;   cfg.rhoFluid = 1.225;      cfg.g = 9.81;
cfg.shapeModel  = 'shape3d';                          % <-- 3D model (twist)
cfg.aero        = struct('C_span', 0.2, 'C_span_torque', 0.7);

% --- Simulation / classification ------------------------------------------
cfg.tspan = [0 6];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();

% --- Animation output ------------------------------------------------------
% Slow motion (playbackSpeed < 1) makes the fast spin watchable; seedZoom tightens
% the follow-cam onto the seed. Videos are written here (kept out of the repo).
cfg.animFps           = 60;
cfg.animPlaybackSpeed = 0.3;
cfg.showSeedVels      = false;   % true -> also draw per-strip wind in the follow-cam
cfg.seedZoom          = 1.0;
cfg.animDir           = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Twist Tests";
if ~exist(cfg.animDir, 'dir'); mkdir(cfg.animDir); end

% --- Base-seed-params + twist amplitudes ----------------------------------
xh = cfg.spanLength/2;   yh = cfg.chordLength/2;   hs = cfg.spanLength/2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

slightTipDeg = 5;    % case 1 tip twist (deg)
largeTipDeg  = 30;   % case 2 tip twist (deg)
asymTipDeg   = 10;   % case 3 single-side tip twist (deg)

%% 1. Slight twist  (anti-symmetric linear, small)
runTwistCase('Twist_slight_5deg', @(z) (deg2rad(slightTipDeg)/hs) .* z, cfg, baseBsp);

%% 2. Large twist  (anti-symmetric linear, large)
runTwistCase('Twist_large_30deg', @(z) (deg2rad(largeTipDeg)/hs) .* z, cfg, baseBsp);

%% 3. Small asymmetric twist  (only the +z half-wing twisted)
runTwistCase('Twist_asym_singleSide_10deg', @(z) (z > 0) .* (deg2rad(asymTipDeg)/hs) .* z, cfg, baseBsp);


% =========================================================================
% LOCAL: build a twisted seed, integrate one drop, and render its animation
% =========================================================================
function runTwistCase(name, twistFn, cfg, baseBsp)
% name    : label (video filename + animation title).
% twistFn : theta(z) handle (rad) -> the spanwise twist profile.

    % --- Build the twisted seed + integrate (from rest, level) via the shared
    % runSingleMode (which routes shape3d -> seed6DOFODE3D through seedRHS).
    bsp = baseBsp;   bsp.twist = twistFn;
    r = runSingleMode(char(name), [0; 0; 0], [1; 0; 0; 0], [0; 0; 0], cfg, bsp);

    % --- Combined mode animation (follow-cam shows the twisted strips) ------
    vfile = fullfile(cfg.animDir, [char(name) '.mp4']);
    animateModeTrajectory(r.t, r.x, struct( ...
        'videoFile',     char(vfile), ...
        'fps',           cfg.animFps, ...
        'playbackSpeed', cfg.animPlaybackSpeed, ...
        'showSeedVels',  cfg.showSeedVels, ...
        'seedZoom',      cfg.seedZoom, ...
        'title',         char(name), ...
        'seedParams',    r.seedParams));
    fprintf('  -> wrote %s\n', vfile);
end
