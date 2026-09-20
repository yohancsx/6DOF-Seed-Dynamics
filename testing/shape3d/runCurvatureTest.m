%% Curvature animation test -- three curved-seed drops, each rendered to an mp4
% Companion to runTwistTest, but for the CURVATURE (dihedral) path of the
% 3D-shape model. Each case drops a CENTRED-CoM seed with a different spanwise
% curvature and renders the combined mode animation (3D trajectory + CoM-in-body
% + zoomed follow-cam, which draws the bent strips) to a video, so you can
% regenerate the curvature animations in one run.
%
% Cases:
%   1. Slight curvature  : symmetric bowl, small tip dihedral -> nearly flat.
%   2. Large curvature   : symmetric bowl, large tip dihedral -> deep bowl; watch
%                          the descent and cone rather than the spin.
%   3. Small asymmetric  : only the +s half-wing bent (single-side) -> the bent
%      curvature           half loads differently from the flat half.
%
% Unlike twist, curvature moves the strip centres OUT of the body x-z plane, so
% the CoM leaves the plane too -- that out-of-plane CoM offset is the signature
% to look for in the CoM-in-body panel.
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
cfg.shapeModel  = 'shape3d';                          % <-- 3D model (curvature)
cfg.aero        = struct('C_span', 0.2);   % no C_span_torque: the 3D model has no span torque
% NOTE: inert here -- setupSeedShape3D switches the planar span-force hack OFF for
% a CURVED seed (its tilted strips supply that force from geometry, so leaving the
% hack on would double-count). See runCurvatureSuite for the flat-vs-0deg caveat.

% --- Simulation / classification ------------------------------------------
cfg.tspan = [0 6];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();

% --- Animation output ------------------------------------------------------
% Curved seeds settle more slowly than twisted ones and spin less, so the default
% playback is less aggressive than runTwistTest's. Videos are kept out of the repo.
cfg.animFps           = 60;
cfg.animPlaybackSpeed = 0.5;
cfg.showSeedVels      = false;   % true -> also draw per-strip wind in the follow-cam
cfg.seedZoom          = 1.0;
cfg.animDir           = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Curvature Tests";
if ~exist(cfg.animDir, 'dir'); mkdir(cfg.animDir); end

% --- Base-seed-params + curvature amplitudes ------------------------------
xh = cfg.spanLength/2;   yh = cfg.chordLength/2;   hs = cfg.spanLength/2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

slightTipDeg = 5;    % case 1 tip dihedral (deg)
largeTipDeg  = 35;   % case 2 tip dihedral (deg)
asymTipDeg   = 15;   % case 3 single-side tip dihedral (deg)

%% 1. Slight curvature  (symmetric bowl, small)
runCurvatureCase('Curve_slight_5deg', @(s) (deg2rad(slightTipDeg)/hs) .* s, cfg, baseBsp);

%% 2. Large curvature  (symmetric bowl, large)
runCurvatureCase('Curve_large_35deg', @(s) (deg2rad(largeTipDeg)/hs) .* s, cfg, baseBsp);

%% 3. Small asymmetric curvature  (only the +s half-wing bent)
runCurvatureCase('Curve_asym_singleSide_15deg', @(s) (s > 0) .* (deg2rad(asymTipDeg)/hs) .* s, cfg, baseBsp);


% =========================================================================
% LOCAL: build a curved seed, integrate one drop, and render its animation
% =========================================================================
function runCurvatureCase(name, curvatureFn, cfg, baseBsp)
% name        : label (video filename + animation title).
% curvatureFn : phi(s) handle (rad) -> the spanwise dihedral tangent-angle profile,
%               where s is ARC LENGTH along the flat span.

    % --- Build the curved seed + integrate (from rest, level) via the shared
    % runSingleMode (which routes shape3d -> seed6DOFODE3D through seedRHS).
    bsp = baseBsp;   bsp.curvature = curvatureFn;
    r = runSingleMode(char(name), [0; 0; 0], [1; 0; 0; 0], [0; 0; 0], cfg, bsp);
    fprintf('  spanForce=%d  CoM_y=%+.3e m  (out-of-plane CoM is curvature''s signature)\n', ...
            r.seedParams.enableSpanForce, r.seedParams.massParams.com_t(2,1));

    % --- Combined mode animation (follow-cam shows the bent strips) ---------
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
