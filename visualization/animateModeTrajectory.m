function animateModeTrajectory(t, x, opts)
% ANIMATEMODETRAJECTORY  Multi-panel, mode-coloured flight animation for ANY run.
%
% A combined figure that makes the flight modes of a single drop legible at a
% glance, and how they change over time:
%
%   LEFT   (large): the 3D CoM trajectory in the inertial frame, the path
%                   coloured by the flight mode the seed is in at each moment,
%                   with a moving dot marking the current position.
%   TOP-RIGHT     : the CoM location WITHIN THE SEED BODY -- how the centre of
%                   mass moves around the body (span horizontal, chord vertical,
%                   the seed's geometric centre at the origin), with the seed
%                   planform outlined and a moving dot. Requires opts.seedParams;
%                   without it, this panel falls back to the inertial bird's-eye
%                   ground track.
%   BOT-RIGHT     : a zoomed follow-cam of the actual seed (its planform + local
%                   wind vectors) rendered at its true orientation each frame, so
%                   you can see it spinning / gliding / tumbling up close. Shown
%                   only when opts.seedParams is supplied.
%
% To actually resolve a fast spin, raise .fps (more/ smoother frames) and drop
% .playbackSpeed below 1 (slow motion) -- see the frame-rate options below.
%
% The mode label per instant comes from classifyModeTimeline (a sliding-window
% application of the same classifier used everywhere else), so the colours match
% the phase map and the rest of the toolchain. Works on any trajectory; nothing
% here is specific to the moving-CoM test.
%
% INPUTS
%   t    : Nx1 time vector (s).
%   x    : Nx13 state history [r(3) q(4) v(3) omega(3)] (ode45 output on
%          seed6DOFODE; rows = time). r/v inertial, omega body.
%   opts : (optional) struct:
%     -- output / video --
%     .videoFile     output path (''/[] -> no video, just show)  (default '')
%     .fps           frame rate of the WHOLE animation: sets both the temporal
%                    sampling (frames = duration*fps) and the base video rate.
%                    Raise for smoother motion, lower for a lighter file (30)
%     .playbackSpeed real-time multiplier on the video rate: 1 = real time,
%                    <1 = slow motion (to watch the spin), >1 = fast (default 1)
%     .figSize       [w h] base figure px  (default [1800 760] w/ seed panels,
%                    else [1280 720])
%     .resolution    capture DPI (scales figSize by res/96)      (default 120)
%     -- content --
%     .title        overall figure title                        (default '')
%     .view         [az el] camera for the 3D panel             (default [40 18])
%     .eventTimes   vector of times to mark (e.g. CoM moves), drawn as ticks on
%                   the trajectory and body-CoM paths                (default [])
%     .eventLabels  cellstr matching eventTimes (optional)       (default {})
%     .lineWidth    trajectory line width                        (default 2)
%     .dotSize      moving-marker size (points^2 for the dots)   (default 90)
%     .modeOpts     struct forwarded to classifyModeTimeline     (default struct())
%     .modeTimeline precomputed classifyModeTimeline output; skips reclassifying
%                   (default [] -> computed here)
%     .trailOnly    true -> paths grow with the dot (past only); false -> full
%                   path shown from the start                    (default false)
%     .seedParams   full seedParams (from setupSeedShapeAndMass). When given, the
%                   top-middle panel shows the CoM's path within the seed body
%                   using massParams.com_t AND the right panel shows the zoomed
%                   seed follow-cam; without it, the top panel shows the inertial
%                   ground track and the right panel is omitted   (default [])
%     .seedZoom     zoomed-view half-window = seedZoom * longest seed dim (2.0)
%     .showSeedVels draw per-strip local wind in the zoomed view  (default true)
%     .shapeOpts    struct forwarded to visualizeSeedShape        (default struct())
%     .velOpts      struct forwarded to visualizeSeedLocalVels    (default struct())
%
% OUTPUT
%   none -- shows the figure and, if .videoFile is set, writes a video.

% =========================================================================
% 0. DEFAULTS
% =========================================================================
if nargin < 3 || isempty(opts); opts = struct(); end
opts = setDefault(opts, 'videoFile',     '');
opts = setDefault(opts, 'fps',           30);    % frame rate of the whole animation
opts = setDefault(opts, 'playbackSpeed', 1);     % >1 faster, <1 slow-motion (see the spin)
opts = setDefault(opts, 'resolution',    120);
opts = setDefault(opts, 'title',         '');
opts = setDefault(opts, 'view',          [40 18]);
opts = setDefault(opts, 'eventTimes',    []);
opts = setDefault(opts, 'eventLabels',   {});
opts = setDefault(opts, 'lineWidth',     2);
opts = setDefault(opts, 'dotSize',       90);
opts = setDefault(opts, 'modeOpts',      struct());
opts = setDefault(opts, 'modeTimeline',  []);
opts = setDefault(opts, 'trailOnly',     false);
opts = setDefault(opts, 'seedParams',    []);    % enables the body-frame CoM + zoomed seed panels
opts = setDefault(opts, 'seedZoom',      1.5);   % zoomed-view half-window = seedZoom * longest seed dim
opts = setDefault(opts, 'showSeedVels',  true);  % draw per-strip local wind in the zoomed view
opts = setDefault(opts, 'shapeOpts',     struct());  % forwarded to visualizeSeedShape
opts = setDefault(opts, 'velOpts',       struct());  % forwarded to visualizeSeedLocalVels

% The body-frame CoM panel and the zoomed seed follow-cam both need the seed
% geometry, so both appear only when seedParams is supplied.
bodyMode = ~isempty(opts.seedParams);
if ~isfield(opts, 'figSize')
    if bodyMode; opts.figSize = [1360 780]; else; opts.figSize = [1180 560]; end
end

t = t(:);
r = x(:, 1:3);   q = x(:, 4:7);   % world CoM position, body->world quaternion

% =========================================================================
% 1. MODE TIMELINE (per-sample label) + palette
% =========================================================================
if isempty(opts.modeTimeline)
    tl = classifyModeTimeline(t, x, opts.modeOpts);
else
    tl = opts.modeTimeline;
end
idx      = tl.idx(:);      % Nx1 mode index per sample
modeList = tl.list;
modeCol  = tl.col;

% =========================================================================
% 2. PLOT-SPACE PATHS
%    3D panel: world Y-up remap [X;Z;Y] (as in the other visualizers).
%    Inertial top view (fallback): horizontal plane = (world X, world Z).
% =========================================================================
P3   = [r(:,1).'; r(:,3).'; r(:,2).'];   % 3xN plot coords for the 3D panel
Ptop = [r(:,1).'; r(:,3).'];             % 2xN (worldX, worldZ) inertial ground track

% Body-frame CoM path (top-right panel), if seedParams was supplied. The CoM in
% body-DATUM coords (massParams.com_t) is recentred on the seed's geometric
% (wing-centroid) so the geometric centre sits at the origin; plotted as
% (spanwise = body z, chordwise = body x). Body z = drawing x, body x = drawing y.
if bodyMode
    sp    = opts.seedParams;
    comB  = sp.massParams.com_t;                 % 3xN_mass, body-datum (m)
    tS    = sp.massParams.tSamples(:);           % N_mass x 1 (s)
    [cSpanCtr, cChordCtr] = centroid(sp.baseSeedParams.seedShape);  % drawX=span, drawY=chord
    comSpan  = comB(3, :).' - cSpanCtr;          % body z minus span-centre
    comChord = comB(1, :).' - cChordCtr;         % body x minus chord-centre
    % Seed planform outline, recentred, in (span, chord)
    V        = sp.baseSeedParams.seedShape.Vertices;   % Kx2 (drawX, drawY), NaN-separated
    outSpan  = V(:, 1) - cSpanCtr;
    outChord = V(:, 2) - cChordCtr;
end

% =========================================================================
% 3. FRAME RESAMPLE onto a uniform fps grid (for smooth constant-rate video)
% =========================================================================
dur = max(t(end) - t(1), eps);
M   = max(2, round(dur * opts.fps));
ti  = linspace(t(1), t(end), M);
P3i   = interp1(t, P3.',   ti).';        % 3xM
Ptopi = interp1(t, Ptop.', ti).';        % 2xM (inertial fallback)
idxi  = interp1(t, idx,    ti, 'nearest', 'extrap');   % 1xM mode per frame
idxi  = idxi(:);                          % column, for consistent indexing

if bodyMode
    if numel(tS) < 2                      % constant CoM -> a single fixed point
        Pbody = repmat([comSpan(1); comChord(1)], 1, M);
    else
        Pbody = [interp1(tS, comSpan,  ti, 'linear', 'extrap'); ...
                 interp1(tS, comChord, ti, 'linear', 'extrap')];   % 2xM (span; chord)
    end

    % --- Extra per-frame data for the zoomed seed follow-cam --------------
    comW_i = interp1(t, r,            ti).';        % 3xM world CoM position
    vW_i   = interp1(t, x(:, 8:10),   ti).';        % 3xM world CoM velocity
    om_i   = interp1(t, x(:, 11:13),  ti).';        % 3xM body angular velocity
    qraw   = interp1(t, q,            ti);          % Mx4 quaternion (linear)
    q_i    = (qraw ./ vecnorm(qraw, 2, 2)).';       % 4xM, renormalised
    if numel(tS) < 2
        comB_i  = repmat(comB(:, 1), 1, M);         % 3xM body CoM (constant)
        massIdx = ones(1, M);
    else
        comB_i  = interp1(tS, comB.', ti).';        % 3xM body CoM
        massIdx = round(interp1(tS, 1:numel(tS), ti, 'nearest', 'extrap'));
        massIdx = min(max(massIdx, 1), numel(tS));  % clamp to valid range
    end
    % Longest planform dimension sets the zoom window.
    [bx, by]   = boundingbox(sp.baseSeedParams.seedShape);
    longestDim = max(bx(2) - bx(1), by(2) - by(1));
end

% =========================================================================
% 4. FIGURE + STATIC CONTENT (built once)
% =========================================================================
resScale  = opts.resolution / 96;
figPx     = round(opts.figSize .* resScale);
fig = figure('Color', 'w', 'Position', [80, 80, figPx(1), figPx(2)]);
% With seed geometry: 2x2 grid -- left 3D (spans both rows) | top-right
% CoM-in-body | bottom-right zoomed seed follow-cam. Without it: 1x2 --
% 3D trajectory | inertial ground track.
if bodyMode
    tlo = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    tlo = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end
if ~isempty(opts.title)
    title(tlo, opts.title, 'FontWeight', 'bold', 'Interpreter', 'none');
end

% --- LEFT: 3D trajectory, coloured by mode ------------------------------
if bodyMode
    ax1 = nexttile(tlo, 1, [2 1]);   % spans both rows of the 2x2
else
    ax1 = nexttile(tlo, 1);          % left of the 1x2
end
hold(ax1, 'on');  grid(ax1, 'on');  box(ax1, 'on');
if ~opts.trailOnly
    drawColoredPath3(ax1, P3, idx, modeCol, opts.lineWidth, '-');
end
% one dummy handle per PRESENT mode -> a clean single-entry legend
present = unique(idx(:)).';
hLeg = gobjects(1, numel(present));
for jj = 1:numel(present)
    hLeg(jj) = plot3(ax1, nan, nan, nan, '-', 'Color', modeCol(present(jj), :), ...
                     'LineWidth', 3, 'DisplayName', modeList{present(jj)});
end
markEvents3(ax1, opts, t, P3);
hGrow3 = plot3(ax1, nan, nan, nan, '-', 'Color', [0.2 0.2 0.2], ...
               'LineWidth', opts.lineWidth, 'HandleVisibility', 'off'); % trailOnly path
hDot3 = scatter3(ax1, nan, nan, nan, opts.dotSize, [0.5 0.5 0.5], 'filled', ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel(ax1, 'World X (m)');  ylabel(ax1, 'World Z (m)');  zlabel(ax1, 'World Y (m) -- up');
view(ax1, opts.view);
setBounds3(ax1, P3);
legend(ax1, hLeg, 'Location', 'northeast', 'Interpreter', 'none');
title(ax1, '3D trajectory (coloured by flight mode)');

% --- TOP-RIGHT: CoM within the seed body (or inertial ground track) ------
ax2 = nexttile(tlo, 2, [1 1]);
hold(ax2, 'on');  grid(ax2, 'on');  box(ax2, 'on');
if bodyMode
    % seed planform outline (recentred on the geometric centre)
    plot(ax2, outSpan, outChord, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(ax2, 0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1, 'HandleVisibility', 'off'); % geom centre
    if ~opts.trailOnly
        drawColoredPath2(ax2, Pbody, idxi, modeCol, 1.5, '--');   % dashed CoM path in body
    end
    markEventsBody(ax2, opts, Pbody, ti);
    hGrow2 = plot(ax2, nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    hDot2  = scatter(ax2, nan, nan, opts.dotSize, [0.5 0.5 0.5], 'filled', ...
                     'MarkerEdgeColor', 'k', 'LineWidth', 1);
    axis(ax2, 'equal');
    xlabel(ax2, 'Spanwise  (body z, m)');  ylabel(ax2, 'Chordwise  (body x, m)');
    title(ax2, 'CoM within seed body (geometric centre at origin)');
    setBodyBounds(ax2, outSpan, outChord, Pbody);
else
    if ~opts.trailOnly
        drawColoredPath2(ax2, Ptop, idx, modeCol, 1.5, '--');   % dashed ground track
    end
    markEvents2(ax2, opts, t, Ptop);
    hGrow2 = plot(ax2, nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    hDot2  = scatter(ax2, nan, nan, opts.dotSize, [0.5 0.5 0.5], 'filled', ...
                     'MarkerEdgeColor', 'k', 'LineWidth', 1);
    axis(ax2, 'equal');
    xlabel(ax2, 'World X (m)');  ylabel(ax2, 'World Z (m)');
    title(ax2, 'Top view (bird''s-eye): CoM ground track');
end

% --- BOT-RIGHT: zoomed seed follow-cam (drawn per frame; static props) ---
if bodyMode
    ax4 = nexttile(tlo, 4, [1 1]);
    hold(ax4, 'on');  grid(ax4, 'on');  box(ax4, 'on');
    daspect(ax4, [1 1 1]);              % true seed proportions
    view(ax4, opts.view);
    xlabel(ax4, 'World X (m)');  ylabel(ax4, 'World Z (m)');  zlabel(ax4, 'World Y (m) -- up');
    title(ax4, 'Seed (zoomed, follows CoM)');
    nTrail = max(2, round(1.5 * opts.fps));   % frames of trailing CoM path shown
end

% =========================================================================
% 5. OPTIONAL VIDEO WRITER
% =========================================================================
writeVid = ~isempty(opts.videoFile);
if writeVid
    [~, ~, ext] = fileparts(opts.videoFile);
    if strcmpi(ext, '.mp4'); profile = 'MPEG-4'; else; profile = 'Motion JPEG AVI'; end
    vw = VideoWriter(opts.videoFile, profile);
    vw.FrameRate = max(0.1, opts.fps * opts.playbackSpeed);   % playbackSpeed<1 -> slow motion
    open(vw);
end

% =========================================================================
% 6. ANIMATION LOOP -- only the dynamic objects change each frame
% =========================================================================
lastPct  = -1;
frameHW  = [];      % locked frame size (even dims) so every writeVideo matches
for f = 1:M
    cm = modeCol(idxi(f), :);

    % current-position dots (colour = current mode)
    set(hDot3, 'XData', P3i(1,f), 'YData', P3i(2,f), 'ZData', P3i(3,f), ...
               'CData', cm);
    if bodyMode
        set(hDot2, 'XData', Pbody(1,f), 'YData', Pbody(2,f), 'CData', cm);
    else
        set(hDot2, 'XData', Ptopi(1,f), 'YData', Ptopi(2,f), 'CData', cm);
    end

    % trailOnly: grow the paths behind the dot
    if opts.trailOnly
        upto = ti <= ti(f);
        set(hGrow3, 'XData', P3i(1,upto), 'YData', P3i(2,upto), 'ZData', P3i(3,upto));
        if bodyMode
            set(hGrow2, 'XData', Pbody(1,upto), 'YData', Pbody(2,upto));
        else
            set(hGrow2, 'XData', Ptopi(1,upto), 'YData', Ptopi(2,upto));
        end
    end

    title(ax1, sprintf('3D trajectory  --  t = %.2f s   |   mode: %s', ...
                       ti(f), modeList{idxi(f)}), 'Interpreter', 'none');

    % --- Zoomed seed follow-cam: full redraw of the seed at this pose -----
    if bodyMode
        Rk       = quatToRotm(q_i(:, f));            % body->world
        comWk    = comW_i(:, f);                     % world CoM (3x1)
        comBk    = comB_i(:, f);                     % body CoM  (3x1)
        posDatum = comWk - Rk * comBk;               % place datum so CoM lands on comWk
        set(fig, 'CurrentAxes', ax4);   cla(ax4);    % redraw this panel only

        % faint recent CoM trail (so translation/glide direction is visible)
        a0 = max(1, f - nTrail);
        plot3(ax4, P3i(1, a0:f), P3i(2, a0:f), P3i(3, a0:f), ...
              '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

        so = opts.shapeOpts;   so.tIndex = massIdx(f);
        visualizeSeedShape(opts.seedParams, Rk, posDatum, so);
        if opts.showSeedVels
            vo = opts.velOpts;   vo.tIndex = massIdx(f);
            slv = computeSeedLocalVel(opts.seedParams, comBk, q_i(:, f), vW_i(:, f), om_i(:, f));
            visualizeSeedLocalVels(opts.seedParams, slv, Rk, posDatum, vo);
        end

        % follow window centred on the seed (plot coords: X, Z, Y-up)
        c = [comWk(1); comWk(3); comWk(2)];   h = opts.seedZoom * longestDim;
        xlim(ax4, [c(1)-h, c(1)+h]);  ylim(ax4, [c(2)-h, c(2)+h]);  zlim(ax4, [c(3)-h, c(3)+h]);
    end

    drawnow;
    if writeVid
        cd = getframe(fig).cdata;
        if isempty(frameHW)
            frameHW = [size(cd,1), size(cd,2)] - mod([size(cd,1), size(cd,2)], 2); % even
        end
        writeVideo(vw, matchFrameSize(cd, frameHW(1), frameHW(2)));
    end

    pct = floor(100 * f / M);
    if pct >= lastPct + 10
        fprintf('animateModeTrajectory: %d%% (%d/%d frames)\n', pct, f, M);
        lastPct = pct;
    end
end

if writeVid
    close(vw);
    fprintf('animateModeTrajectory: wrote %d frames to "%s".\n', M, opts.videoFile);
end
end   % animateModeTrajectory


% =========================================================================
% LOCAL: draw a 3xN path as mode-coloured runs (one line per contiguous mode)
% =========================================================================
function drawColoredPath3(ax, P, idx, col, lw, style)
    idx = idx(:);
    n  = size(P, 2);
    ch = [1; find(diff(idx) ~= 0) + 1; n + 1];   % run-start indices + sentinel
    for r = 1:numel(ch) - 1
        a  = ch(r);   b = ch(r+1) - 1;
        aa = max(a - 1, 1);                       % overlap one point to connect runs
        plot3(ax, P(1, aa:b), P(2, aa:b), P(3, aa:b), style, ...
              'Color', col(idx(a), :), 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end


% =========================================================================
% LOCAL: draw a 2xN path as mode-coloured runs
% =========================================================================
function drawColoredPath2(ax, P, idx, col, lw, style)
    idx = idx(:);
    n  = size(P, 2);
    ch = [1; find(diff(idx) ~= 0) + 1; n + 1];
    for r = 1:numel(ch) - 1
        a  = ch(r);   b = ch(r+1) - 1;
        aa = max(a - 1, 1);
        plot(ax, P(1, aa:b), P(2, aa:b), style, ...
             'Color', col(idx(a), :), 'LineWidth', lw, 'HandleVisibility', 'off');
    end
end


% =========================================================================
% LOCAL: event markers (ticks on the trajectory / body-CoM paths)
% =========================================================================
function markEvents3(ax, opts, t, P3)
    if isempty(opts.eventTimes); return; end
    for e = opts.eventTimes(:).'
        [~, k] = min(abs(t - e));
        plot3(ax, P3(1,k), P3(2,k), P3(3,k), 'kx', 'MarkerSize', 9, ...
              'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

function markEvents2(ax, opts, t, Ptop)
    if isempty(opts.eventTimes); return; end
    for e = opts.eventTimes(:).'
        [~, k] = min(abs(t - e));
        plot(ax, Ptop(1,k), Ptop(2,k), 'kx', 'MarkerSize', 9, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

function markEventsBody(ax, opts, Pbody, ti)
% Mark each CoM-move instant at its body-frame position (span, chord).
    if isempty(opts.eventTimes); return; end
    for e = opts.eventTimes(:).'
        [~, k] = min(abs(ti - e));
        plot(ax, Pbody(1,k), Pbody(2,k), 'kx', 'MarkerSize', 9, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end


% =========================================================================
% LOCAL: body-frame panel bounds -- enclose the seed outline AND the CoM
% excursion (which can run outside the body), padded; equal aspect is set by
% the caller so the seed's true wide/thin proportions show
% =========================================================================
function setBodyBounds(ax, outSpan, outChord, Pbody)
    xs = [outSpan(:); Pbody(1, :).'];   xs = xs(isfinite(xs));
    ys = [outChord(:); Pbody(2, :).'];  ys = ys(isfinite(ys));
    padx = 0.10 * (max(xs) - min(xs)) + eps;
    pady = 0.10 * (max(ys) - min(ys)) + eps;
    xlim(ax, [min(xs) - padx, max(xs) + padx]);
    ylim(ax, [min(ys) - pady, max(ys) + pady]);
end


% =========================================================================
% LOCAL: 3D axis bounds enclosing the whole path (padded), no equal aspect
% (a falling seed drops far more than it drifts; equal aspect would squash the
% horizontal structure into a line, so the box is filled to reveal the mode shape)
% =========================================================================
function setBounds3(ax, P3)
    mn = min(P3, [], 2);   mx = max(P3, [], 2);
    pad = 0.05 * max(mx - mn, [], 'all') + eps;
    span = max(mx - mn, 1e-3);
    xlim(ax, [mn(1)-0.1*span(1)-pad, mx(1)+0.1*span(1)+pad]);
    ylim(ax, [mn(2)-0.1*span(2)-pad, mx(2)+0.1*span(2)+pad]);
    zlim(ax, [mn(3)-0.05*span(3)-pad, mx(3)+0.05*span(3)+pad]);
end


% =========================================================================
% LOCAL: crop/pad a captured frame to a fixed HxW (guards against getframe
% returning off-by-a-pixel sizes when the 3D follow-cam re-renders, which
% would make writeVideo reject the frame)
% =========================================================================
function cd = matchFrameSize(cd, H, W)
    h = size(cd, 1);
    if h > H
        cd = cd(1:H, :, :);
    elseif h < H
        cd = cat(1, cd, 255 * ones(H - h, size(cd, 2), 3, 'uint8'));
    end
    w = size(cd, 2);
    if w > W
        cd = cd(:, 1:W, :);
    elseif w < W
        cd = cat(2, cd, 255 * ones(H, W - w, 3, 'uint8'));
    end
end


% =========================================================================
% LOCAL: set default struct field if absent
% =========================================================================
function s = setDefault(s, field, value)
    if ~isfield(s, field); s.(field) = value; end
end
