function animateSeed(seedParams, traj, opts)
% ANIMATESEED  Render a seed trajectory to a video, frame by frame.
%
% Wraps visualizeSeedShape and visualizeSeedLocalVels in a per-frame loop.
% Each frame is drawn and streamed straight into the VideoWriter as it's
% generated -- no per-frame PNG round-trip to disk unless opts.saveFrameImages
% is set, so the common case (video only) skips a redundant write-then-re-read
% of every frame. Designed to be extended: new per-frame overlays (forces,
% torques, ...) slot into renderFrame behind their own opts toggle without
% touching the loop.
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
%     .videoFile       output video path         (default: 'seedAnimation.mp4')
%     .fps             video frame rate           (default: 30)
%     .resolution      output DPI. Frames are captured with getframe at the
%                      figure's actual pixel size, so this scales .figSize by
%                      resolution/96 (96 DPI = native figSize, no scaling) to
%                      approximate what a DPI-based export would have produced
%                                                        (default: 150)
%     .figSize         [w h] BASE figure px, before .resolution scaling
%                                            (default: [900 700])
%     .frameStep       render every Nth STATE as a video frame (default: 1)
%     .saveFrameImages also save individual rendered frames as PNGs, in
%                      addition to the video                    (default: false)
%     .outputFolder    folder for saved PNG frames -- only used when
%                      .saveFrameImages is true       (default: 'anim_frames')
%     .frameImageStep  of the RENDERED frames, save every Nth one as a PNG
%                      when .saveFrameImages is true   (default: 10)
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
%     .fontScale       multiplier on the default axes font size, applied to
%                      the title, tick labels, and axis labels together
%                      (via ax.FontSize)                (default: 1.5)

% =========================================================================
% 0. DEFAULT OPTIONS
% =========================================================================
if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = setDefault(opts, 'videoFile',       'seedAnimation.mp4');
opts = setDefault(opts, 'fps',             30);
opts = setDefault(opts, 'resolution',      150);
opts = setDefault(opts, 'figSize',         [900 700]);
opts = setDefault(opts, 'frameStep',       1);
opts = setDefault(opts, 'saveFrameImages', false);
opts = setDefault(opts, 'outputFolder',    'anim_frames');
opts = setDefault(opts, 'frameImageStep',  10);

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
opts = setDefault(opts, 'fontScale',  1.5);

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

% Fresh output folder, only needed if individual frame images are saved.
if opts.saveFrameImages && ~exist(opts.outputFolder, 'dir')
    mkdir(opts.outputFolder);
end

% =========================================================================
% 2. RENDER LOOP -- stream frames straight into the video writer
% =========================================================================
% One reusable figure AND axes, both created once. Contents are cleared each
% frame with cla(ax) rather than rebuilding the axes from scratch (clf +
% axes()) -- axes creation is one of the more expensive per-frame costs, and
% is entirely avoidable since the axes' static properties (below) don't
% change frame to frame. Kept VISIBLE -- an invisible ('Visible','off')
% figure combined with the transparent/OpenGL-rendered seed patches produced
% black, badly-slow exports; a visible figure renders correctly and lets you
% watch it live.
resScale  = opts.resolution / 96;               % DPI relative to a 96 DPI screen
figPxSize = round(opts.figSize .* resScale);    % actual on-screen/captured size

fig = figure('Color', 'w', ...
             'Position', [100, 100, figPxSize(1), figPxSize(2)]);
ax  = axes('Parent', fig);
hold(ax, 'on');

% Static axes properties -- set ONCE. cla(ax) inside the loop clears plotted
% content (patches/lines/markers) but does not touch axes properties like
% these, or the Title/XLabel/YLabel/ZLabel objects, so repeating them every
% frame (as the previous version did) was pure waste.
daspect(ax, [1 1 1]);                 % equal aspect (cube)
grid(ax, 'on');
view(ax, opts.view);
xlabel(ax, 'World X (m)');
ylabel(ax, 'World Z (m)');
zlabel(ax, 'World Y (m) -- up');
ax.FontSize = get(groot, 'defaultAxesFontSize') * opts.fontScale;   % scales
    % ticks directly; title/axis labels scale with it via their default
    % TitleFontSizeMultiplier/LabelFontSizeMultiplier.

% 'fixed' axis limits don't change per frame either -- set once here.
% 'follow' limits track the seed's position, so those are still set inside
% the loop below.
if strcmpi(opts.axisMode, 'fixed')
    xlim(ax, fixedLimits(1, :)); ylim(ax, fixedLimits(2, :)); zlim(ax, fixedLimits(3, :));
end

% Pick a VideoWriter profile from the file extension. MPEG-4 needs platform
% support (Windows/macOS); fall back to Motion JPEG AVI elsewhere. Opened
% before the loop so each frame can be written as soon as it's captured,
% instead of round-tripping every frame through a PNG on disk first.
[~, ~, ext] = fileparts(opts.videoFile);
if strcmpi(ext, '.mp4')
    profile = 'MPEG-4';
else
    profile = 'Motion JPEG AVI';
end
vw = VideoWriter(opts.videoFile, profile);
vw.FrameRate = opts.fps;
open(vw);

numRenderFrames = numel(frameIdx);
savedFrameFiles = strings(1, 0);   % PNGs actually saved, if any
lastPctPrinted  = -1;              % forces a print at the first 0%+ boundary

for f = 1 : numRenderFrames
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
    cla(ax);   % clear plotted content only -- axes object/properties persist

    renderFrame(seedParams, traj, k, R_b2w, posDatum, comB, q, ...
                trajPlot, toPlot, canDoVels, opts);

    % --- Axis limits: only 'follow' mode changes per frame; 'fixed' was
    % already set once above. -----------------------------------------------
    if strcmpi(opts.axisMode, 'follow')
        lims = computeFollowLimits(toPlot(comW), longestDim, opts.zoomFactor);
        xlim(ax, lims(1, :)); ylim(ax, lims(2, :)); zlim(ax, lims(3, :));
    end
    if isfield(traj, 't')
        title(ax, sprintf('t = %.3f s   (frame %d/%d)', traj.t(k), f, numRenderFrames));
    else
        title(ax, sprintf('frame %d/%d', f, numRenderFrames));
    end

    % --- Write directly to the video (no PNG round-trip) -----------------
    % getframe grabs the already-rendered pixel buffer directly, instead of
    % exportgraphics's print-style re-render -- much faster in a per-frame
    % loop, and doesn't hit the invisible-figure black-render bug.
    drawnow;
    frame = getframe(fig);
    writeVideo(vw, frame.cdata);

    % --- Optionally also save this frame as a PNG -------------------------
    if opts.saveFrameImages && mod(f - 1, opts.frameImageStep) == 0
        fname = fullfile(opts.outputFolder, sprintf('frame_%05d.png', f));
        imwrite(frame.cdata, fname);
        savedFrameFiles(end+1) = fname; %#ok<AGROW> -- opt-in path, rarely grows
    end

    % --- Progress, printed roughly every 10% ------------------------------
    pct = floor(100 * f / numRenderFrames);
    if pct >= lastPctPrinted + 10
        fprintf('animateSeed: %d%% done (%d/%d frames)\n', pct, f, numRenderFrames);
        lastPctPrinted = pct;
    end
end

close(fig);
close(vw);

if opts.saveFrameImages
    fprintf('animateSeed: wrote %d frames to video "%s" (%d frame images saved to "%s").\n', ...
            numRenderFrames, opts.videoFile, numel(savedFrameFiles), opts.outputFolder);
else
    fprintf('animateSeed: wrote %d frames to video "%s".\n', numRenderFrames, opts.videoFile);
end

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