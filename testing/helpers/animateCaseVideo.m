function ok = animateCaseVideo(t, x, sp, name, outDir, aopt)
% ANIMATECASEVIDEO  Write one mode-coloured .mp4 (+ a final-frame .png) for a
%   single case, with the windowing / housekeeping a sweep runner needs.
%
% A thin wrapper around animateModeTrajectory that does the four things every
% caller in a sweep otherwise has to repeat: create the output folder (the
% animator does not), trim the flight to the interesting window, close the
% figure afterwards, and -- the point of it -- NEVER throw. A rendering failure
% on one case must not abort a sweep that has already spent minutes integrating,
% so failures are caught, reported, and returned as a false status.
%
% WHY A WINDOW: animation cost is linear in frames, and a 12 s drop at 30 fps is
% 360 rendered frames PER CASE. The settled behaviour -- the thing being compared
% across a sweep -- lives at the end of the flight, so by default only the last
% aopt.window seconds are drawn. Set window = Inf for the whole flight.
%
% INPUTS
%   t, x   : ode45 outputs, Nx1 time and Nx13 state [r q v omega].
%   sp     : the seedParams that were integrated; enables the body-CoM panel and
%            the zoomed seed follow-cam. Pass [] for the trajectory-only layout.
%   name   : case name -- the video/figure title and the file stem.
%   outDir : folder to write into (created if absent).
%   aopt   : (optional) struct:
%     .window        seconds of flight to animate, counted back from the END
%                    (default 4; Inf = whole flight)
%     .fps           frame rate                                   (default 20)
%     .playbackSpeed <1 slows it down enough to see a fast spin    (default 0.25)
%     .showSeedVels  per-strip local wind arrows in the zoom panel (default false)
%     .savePng       also export the final frame as <name>.png     (default true)
%     .eventTimes    times to mark on the paths (e.g. CoM moves)   (default [])
%     .blind         no mode colours/labels (for by-eye labelling)  (default false)
%     .topPanel      'auto' | 'body' | 'groundTrack'               (default 'auto')
%
% OUTPUT
%   ok : true if the video was written, false if it failed (reason printed).

    if nargin < 6 || isempty(aopt); aopt = struct(); end
    if ~isfield(aopt,'window');        aopt.window        = 4;     end
    if ~isfield(aopt,'fps');           aopt.fps           = 20;    end
    if ~isfield(aopt,'playbackSpeed'); aopt.playbackSpeed = 0.25;  end
    if ~isfield(aopt,'showSeedVels');  aopt.showSeedVels  = false; end
    if ~isfield(aopt,'savePng');       aopt.savePng       = true;  end
    if ~isfield(aopt,'eventTimes');    aopt.eventTimes    = [];    end
    if ~isfield(aopt,'blind');         aopt.blind         = false; end
    if ~isfield(aopt,'topPanel');      aopt.topPanel      = 'auto'; end

    ok = false;
    if isempty(t) || isempty(x) || any(~isfinite(x(:)))
        fprintf(2, '    animation skipped (%s): no finite trajectory\n', name);
        return
    end
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    % --- Trim to the window, keeping at least a handful of samples ----------
    sel = true(size(t));
    if isfinite(aopt.window) && aopt.window > 0
        sel = t >= (t(end) - aopt.window);
        if nnz(sel) < 8; sel = true(size(t)); end
    end
    tw = t(sel);   xw = x(sel, :);

    stem  = matlab.lang.makeValidName(char(name));
    vfile = fullfile(outDir, [stem '.mp4']);
    try
        animateModeTrajectory(tw, xw, struct( ...
            'videoFile',     vfile, ...
            'fps',           aopt.fps, ...
            'playbackSpeed', aopt.playbackSpeed, ...
            'showSeedVels',  aopt.showSeedVels, ...
            'eventTimes',    aopt.eventTimes, ...
            'blind',         aopt.blind, ...
            'topPanel',      aopt.topPanel, ...
            'title',         char(name), ...
            'seedParams',    sp));
        if aopt.savePng
            exportgraphics(gcf, fullfile(outDir, [stem '.png']), 'Resolution', 120);
        end
        close(gcf);
        ok = true;
    catch ME
        fprintf(2, '    animation FAILED (%s): %s\n', name, ME.message);
        if ~isempty(findobj('Type','figure')); close(gcf); end
    end
end
