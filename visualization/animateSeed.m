function animateSeed(seedParams, traj, opts)
% ANIMATESEED  Render a seed trajectory to a video, frame by frame.
%
% Wraps visualizeSeedShape and visualizeSeedLocalVels in a per-frame loop.
% Each frame is drawn, saved as a PNG into a folder, and afterwards all frames
% are composited into a single video file. Designed to be extended: new
% per-frame overlays (forces, torques, etc.) slot into renderFrame behind their
% own opts toggle without touching the loop or compositor.
%
% INPUTS
%   seedParams : struct from setupSeedShapeAndMass.
%   traj       : trajectory struct with per-frame columns (N frames):
%     .comPos_world   - 3xN CoM position in the inertial frame (m)   [required]
%     .quat           - 4xN body->world quaternions [q0;q1;q2;q3]    [required]
%     .comPos_body    - 3xN CoM location in body-datum coords (m)    [optional;
%                       defaults to massParams.com_t, or its first column]
%     .comVel_world   - 3xN CoM velocity, inertial frame (m/s)       [optional;
%                       required only if velocities are drawn]
%     .comOmega_body  - 3xN body-frame angular velocity (rad/s)      [optional;
%                       required only if velocities are drawn]
%   opts       : (optional) struct of settings. Defaults below.
%     -- output / video --
%     .outputFolder    folder for PNG frames     (default: 'anim_frames')
%     .videoFile       output video path         (default: 'seedAnimation.mp4')
%     .fps             video frame rate           (default: 30)
%     .resolution      export DPI                 (default: 150)
%     .figSize         [w h] figure px            (default: [900 700])
%     .frameStep       render every Nth frame     (default: 1)
%     .keepFrames      keep PNGs after compositing (default: true)
%     .visibleFigure   show the figure while rendering (default: false)
%     -- content toggles --
%     .showShape       draw the seed each frame   (default: true)
%     .showVels        draw local velocities      (default: true)
%     .showTrajTotal   draw full CoM path (faint) (default: true)
%     .showTrajPartial draw CoM path up to frame  (default: true)
%     .shapeOpts       struct forwarded to visualizeSeedShape     (default: struct())
%     .velOpts         struct forwarded to visualizeSeedLocalVels (default: struct())
%     -- trajectory style --
%     .colorTrajTotal   full-path colour   (default: [0.7 0.7 0.7] light grey)
%     .colorTrajPartial trail colour       (default: [0.85 0.10 0.10] red)
%     .trajLineWidth                       (default: 1.5)
%     -- axes --
%     .axisMode        'follow' (zoom on seed) or 'fixed' (whole path)
%                                            (default: 'follow')
%     .zoomFactor      half-window = zoomFactor * longest seed dimension,
%                      so the seed fills the view  (default: 2)
%     .view            [az el] camera angle  (default: [35 20])

% =========================================================================
% 0. DEFAULT OPTIONS
% =========================================================================
if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = setDefault(opts, 'outputFolder',    'anim_frames');
opts = setDefault(opts, 'videoFile',       'seedAnimation.mp4');
opts = setDefault(opts, 'fps',             30);
opts = setDefault(opts, 'resolution',      150);
opts = setDefault(opts, 'figSize',         [900 700]);
opts = setDefault(opts, 'frameStep',       1);
opts = setDefault(opts, 'keepFrames',      true);
opts = setDefault(opts, 'visibleFigure',   false);

opts = setDefault(opts, 'showShape',       true);
opts = setDefault(opts, 'showVels',        true);
opts = setDefault(opts, 'showTrajTotal',   true);
opts = setDefault(opts, 'showTrajPartial', true);
opts = setDefault(opts, 'shapeOpts',       struct());
opts = setDefault(opts, 'velOpts',         struct());

opts = setDefault(opts, 'colorTrajTotal',   [0.70 0.70 0.70]);
opts = setDefault(opts, 'colorTrajPartial', [0.85 0.10 0.10]);
opts = setDefault(opts, 'trajLineWidth',    1.5);

opts = setDefault(opts, 'axisMode',  'follow');
opts = setDefault(opts, 'zoomFactor', 2);
opts = setDefault(opts, 'view',      [35 20]);

% =========================================================================
% 1. PRELIMINARIES
% =========================================================================
numFrames = size(traj.comPos_world, 2);
frameIdx  = 1 : opts.frameStep : numFrames;   % which states become frames

% World->plot mapping (Y-up), identical to the visualizers.
toPlot = @(w) [w(1,:); w(3,:); w(2,:)];

% Longest seed dimension (for axis zoom). Taken from the source polyshape:
% drawing x = body z (span), drawing y = body x (chord).
[xv, yv]   = boundingbox(seedParams.baseSeedParams.seedShape);
longestDim = max(xv(2) - xv(1), yv(2) - yv(1));

% Full CoM path in plot coordinates (used for both trajectory and fixed axes).
trajPlot = toPlot(traj.comPos_world);   % 3xN

% Whether velocity drawing is possible/requested.
canDoVels = opts.showVels ...
            && isfield(traj, 'comVel_world') ...
            && isfield(traj, 'comOmega_body');
if opts.showVels && ~canDoVels
    warning('animateSeed:noVelData', ...
        'showVels is true but traj lacks comVel_world/comOmega_body; skipping velocities.');
end

% Precompute fixed axis limits once (only used in 'fixed' mode).
fixedLimits = computeFixedLimits(trajPlot, longestDim);

% Fresh output folder.
if ~exist(opts.outputFolder, 'dir')
    mkdir(opts.outputFolder);
end

% =========================================================================
% 2. RENDER LOOP -- one PNG per frame
% =========================================================================
% One reusable figure (cleared each frame) to avoid leaking handles.
fig = figure('Color', 'w', ...
             'Position', [100, 100, opts.figSize(1), opts.figSize(2)]);
if ~opts.visibleFigure
    set(fig, 'Visible', 'off');
end

frameFiles = strings(1, numel(frameIdx));   % saved filenames, in order

for f = 1 : numel(frameIdx)
    k = frameIdx(f);   % index into the trajectory arrays

    % --- Per-frame pose -------------------------------------------------
    q          = traj.quat(:, k);
    R_b2w      = quatToRotm(q);                 % body->world
    comW       = traj.comPos_world(:, k);       % CoM in world (3x1)
    comB       = getComBody(traj, seedParams, k, numFrames);   % CoM in body
    
    % The visualizers position the body-DATUM origin at 'pos'. We want the
    % CoM to land on comW, so place the datum origin accordingly:
    %   comW = R_b2w * comB + posDatum  ->  posDatum = comW - R_b2w*comB
    posDatum = comW - R_b2w * comB; %could be a source of bugs

    % --- Draw the frame -------------------------------------------------
    clf(fig);
    ax = axes('Parent', fig);
    hold(ax, 'on');

    renderFrame(seedParams, traj, k, R_b2w, posDatum, comB, q, ...
                trajPlot, toPlot, canDoVels, opts);

    % --- Axes: limits, aspect, view, labels -----------------------------
    if strcmpi(opts.axisMode, 'follow')
        lims = computeFollowLimits(toPlot(comW), longestDim, opts.zoomFactor);
    else
        lims = fixedLimits;
    end
    xlim(ax, lims(1, :)); ylim(ax, lims(2, :)); zlim(ax, lims(3, :));
    daspect(ax, [1 1 1]);                 % equal aspect (cube)
    grid(ax, 'on');
    view(ax, opts.view);
    xlabel(ax, 'World X (m)');
    ylabel(ax, 'World Z (m)');
    zlabel(ax, 'World Y (m) -- up');
    if isfield(traj, 't')
        title(ax, sprintf('t = %.3f s   (frame %d/%d)', traj.t(k), f, numel(frameIdx)));
    else
        title(ax, sprintf('frame %d/%d', f, numel(frameIdx)));
    end

    % --- Save frame image ----------------------------------------------
    fname = fullfile(opts.outputFolder, sprintf('frame_%05d.png', f));
    exportgraphics(fig, fname, 'Resolution', opts.resolution);
    frameFiles(f) = fname;
end

close(fig);

% =========================================================================
% 3. COMPOSITE FRAMES INTO A VIDEO
% =========================================================================
% Pick a VideoWriter profile from the file extension. MPEG-4 needs platform
% support (Windows/macOS); fall back to Motion JPEG AVI elsewhere.
[~, ~, ext] = fileparts(opts.videoFile);
if strcmpi(ext, '.mp4')
    profile = 'MPEG-4';
else
    profile = 'Motion JPEG AVI';
end

vw = VideoWriter(opts.videoFile, profile);
vw.FrameRate = opts.fps;
open(vw);
for f = 1 : numel(frameFiles)
    img = imread(frameFiles(f));
    writeVideo(vw, img);
end
close(vw);
% 
% % =========================================================================
% % 4. OPTIONAL CLEANUP
% % =========================================================================
% if ~opts.keepFrames
%     for f = 1 : numel(frameFiles)
%         delete(frameFiles(f));
%     end
% end

fprintf('animateSeed: wrote %d frames to "%s" and video "%s".\n', ...
        numel(frameFiles), opts.outputFolder, opts.videoFile);

end   % animateSeed


% =========================================================================
% LOCAL: render a single frame's contents (the extension point)
% =========================================================================
function renderFrame(seedParams, traj, k, R_b2w, posDatum, comB, q, ...
                     trajPlot, toPlot, canDoVels, opts)
% Draws everything for frame k onto the current axes. Add future overlays
% (forces, torques, ...) here behind their own opts toggle.

    % --- Trajectories ---------------------------------------------------
    if opts.showTrajTotal
        plot3(trajPlot(1, :), trajPlot(2, :), trajPlot(3, :), ...
              '-', 'Color', opts.colorTrajTotal, 'LineWidth', opts.trajLineWidth);
    end
    if opts.showTrajPartial
        plot3(trajPlot(1, 1:k), trajPlot(2, 1:k), trajPlot(3, 1:k), ...
              '-', 'Color', opts.colorTrajPartial, 'LineWidth', opts.trajLineWidth);
    end

    % --- Seed shape -----------------------------------------------------
    if opts.showShape
        visualizeSeedShape(seedParams, R_b2w, posDatum, opts.shapeOpts);
    end

    % --- Local velocities ----------------------------------------------
    if canDoVels
        seedLocalVels = computeSeedLocalVel( ...
            seedParams, comB, q, ...
            traj.comVel_world(:, k), traj.comOmega_body(:, k));
        visualizeSeedLocalVels(seedParams, seedLocalVels, R_b2w, posDatum, opts.velOpts);
    end

    % --- [EXTENSION POINT] future overlays, e.g.:
    %   if opts.showForces
    %       visualizeSeedForces(seedParams, seedForces, R_b2w, posDatum, opts.forceOpts);
    %   end
end


% =========================================================================
% LOCAL: CoM in body-datum coords for frame k (with sensible fallbacks)
% =========================================================================
function comB = getComBody(traj, seedParams, k, numFrames)
    if isfield(traj, 'comPos_body')
        comB = traj.comPos_body(:, k);
    else
        com_t = seedParams.massParams.com_t;
        if size(com_t, 2) == numFrames
            comB = com_t(:, k);        % per-frame body CoM available
        else
            comB = com_t(:, 1);        % constant body CoM
        end
    end
end


% =========================================================================
% LOCAL: follow-mode axis limits (cube centred on the seed)
% =========================================================================
function lims = computeFollowLimits(centerPlot, longestDim, zoomFactor)
% centerPlot : 3x1 seed centre in plot coords. Half-window scales with the
% longest seed dimension so the seed always fills a consistent fraction of view.
    h = zoomFactor * longestDim;
    lims = [centerPlot - h, centerPlot + h];   % 3x2: [min max] per plot axis
end


% =========================================================================
% LOCAL: fixed-mode axis limits (cube enclosing the whole path)
% =========================================================================
function lims = computeFixedLimits(trajPlot, longestDim)
    mins = min(trajPlot, [], 2) - longestDim;   % pad by one seed length
    maxs = max(trajPlot, [], 2) + longestDim;
    center = (mins + maxs) / 2;
    halfRange = max((maxs - mins) / 2);          % largest half-extent -> cube
    lims = [center - halfRange, center + halfRange];   % 3x2
end


% =========================================================================
% LOCAL: set default struct field if absent
% =========================================================================
function s = setDefault(s, field, value)
    if ~isfield(s, field)
        s.(field) = value;
    end
end