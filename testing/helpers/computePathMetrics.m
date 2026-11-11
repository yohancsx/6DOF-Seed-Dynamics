function pm = computePathMetrics(t, pos, opts)
% COMPUTEPATHMETRICS  Position-only descent metrics, for simulated AND measured
%   drops alike.
%
% The tracked drops have one point per frame and no attitude, so the comparison
% with the model uses only what a camera sees. The same function, with the same
% smoothing and the same window in height fallen, is applied to the simulated CoM
% path (resampled at the camera rate by the caller) and to the track. Nothing here
% is model-specific.
%
% THE WINDOW is a band of HEIGHT FALLEN, zDrop = z(1) - z, so both sides are
% measured over the same part of the fall: past the release transient, and down
% to where the shorter of the two paths ends (opts.window).
%
% METRICS (all over the window)
%   descentSpeed  slope of a straight-line fit of zDrop(t) (m/s). The raw points
%                 are fitted, so tracking noise averages out.
%   horizSpeed    mean horizontal speed of the SMOOTHED path (m/s).
%   driftSpeed    net horizontal speed: slopes of straight-line fits of x(t), y(t).
%   glideDeg      atan(driftSpeed / descentSpeed): the mean glide angle from the
%                 vertical (deg). 0 = straight down.
%   helix         the RAW horizontal path minus its linear drift, as a complex
%                 signal w(t) = rx + i ry, fitted with one rotating vector
%                 A exp(i 2 pi f t) by scanning f (both signs) for the largest
%                 |sum w exp(-i 2 pi f t)|:
%                   .radius = |A| (m), .turns = f * duration (signed; + =
%                   anticlockwise seen from above), .period = 1/|f| (s),
%                   .explained (fraction of the residual's variance the rotation
%                   accounts for), .rms (m, smoothed residual excursion), .valid
%                 valid = at least one full turn AND explained >= opts.minExplained.
%                 Tracking noise (~4.5 cm on the real drops) spreads over all
%                 frequencies and explains almost nothing, so it cannot pass -- a
%                 spectral-peak test on the SMOOTHED path did, which is why this one
%                 works on the raw points.
%   tortuosity    horizontal path length / net horizontal displacement (smoothed).
%   transientDrop height fallen before the smoothed descent speed first reaches 90%
%                 of descentSpeed (m), over the whole path.
%   pathClass     in this order: 'glide' (glideDeg >= opts.glideDeg; a gliding seed
%                 may also curve, which its helix numbers still show) | 'helix'
%                 (valid helix) |
%                 'straight' (smoothed residual rms < opts.straightRms, i.e. at
%                 the tracking-noise level) | 'irregular'.
%                 A summary label; the numbers above are the comparison.
%
% INPUTS
%   t    : Nx1 times (s), uniformly sampled.
%   pos  : Nx3 positions (m), z vertically UP.
%   opts : (optional) .window [d0 d1] in zDrop (m; default [1.5, max fallen]),
%          .smoothTime (s; Savitzky-Golay span, default 0.25), .minExplained (0.5),
%          .maxHelixHz (10), .glideDeg (20), .straightRms (0.05 m).
% OUTPUT
%   pm   : struct of the metrics above, plus .window, .nWindow, .duration, .valid
%          (false if the window holds fewer than 20 samples; metrics NaN then).

    if nargin < 3; opts = struct(); end
    t = t(:);   zDrop = pos(1,3) - pos(:,3);
    if ~isfield(opts, 'window') || isempty(opts.window); opts.window = [1.5, max(zDrop)]; end
    if ~isfield(opts, 'smoothTime');   opts.smoothTime   = 0.25; end
    if ~isfield(opts, 'minExplained'); opts.minExplained = 0.5;  end
    if ~isfield(opts, 'maxHelixHz');   opts.maxHelixHz   = 10;   end
    if ~isfield(opts, 'glideDeg');     opts.glideDeg     = 20;   end
    if ~isfield(opts, 'straightRms');  opts.straightRms  = 0.05; end

    nanH = struct('radius', NaN, 'turns', NaN, 'period', NaN, 'rms', NaN, 'explained', NaN, 'valid', false);
    pm = struct('descentSpeed', NaN, 'horizSpeed', NaN, 'driftSpeed', NaN, 'glideDeg', NaN, ...
                'helix', nanH, 'tortuosity', NaN, 'transientDrop', NaN, 'pathClass', "none", ...
                'window', opts.window, 'nWindow', 0, 'duration', NaN, 'valid', false);
    win = zDrop >= opts.window(1) & zDrop <= opts.window(2);
    pm.nWindow = nnz(win);
    if pm.nWindow < 20; return; end

    dt = median(diff(t));
    n  = max(5, 2*floor(opts.smoothTime/dt/2) + 1);         % odd span, >= 5 samples
    ps = smoothdata(pos, 1, 'sgolay', min(n, 2*floor(size(pos,1)/2) - 1));
    vel = [gradient(ps(:,1), dt), gradient(ps(:,2), dt), gradient(ps(:,3), dt)];

    tw = t(win);   tw0 = tw - tw(1);
    pm.duration = tw(end) - tw(1);
    c = polyfit(tw0, zDrop(win), 1);              pm.descentSpeed = c(1);
    pm.horizSpeed = mean(hypot(vel(win,1), vel(win,2)));
    cx = polyfit(tw0, pos(win,1), 1);   cy = polyfit(tw0, pos(win,2), 1);
    pm.driftSpeed = hypot(cx(1), cy(1));
    pm.glideDeg   = atand(pm.driftSpeed / max(pm.descentSpeed, eps));

    % --- helix: one rotating vector fitted to the drift-free horizontal path --
    h  = nanH;
    sx = ps(win,1) - polyval(cx, tw0);   sy = ps(win,2) - polyval(cy, tw0);
    h.rms = sqrt(mean(sx.^2 + sy.^2));                    % smoothed excursion
    wz = (pos(win,1) - polyval(cx, tw0)) + 1i*(pos(win,2) - polyval(cy, tw0));
    wz = wz - mean(wz);
    fMax = min(opts.maxHelixHz, 0.25/dt);
    fg   = (1/pm.duration : 1/(8*pm.duration) : fMax).';  % >= one turn in the window
    if ~isempty(fg)
        fg = [fg; -fg];
        E  = exp(-1i*2*pi*tw0(:).' .* fg);                % nf x nw
        A  = (E * wz) / numel(wz);                        % rotating-vector amplitudes
        [~, j] = max(abs(A));
        h.radius = abs(A(j));   h.turns = fg(j) * pm.duration;   h.period = 1/abs(fg(j));
        h.explained = numel(wz) * abs(A(j))^2 / max(sum(abs(wz).^2), eps);
        h.valid = abs(h.turns) >= 1 && h.explained >= opts.minExplained;
    end
    pm.helix = h;

    seg = hypot(diff(ps(win,1)), diff(ps(win,2)));
    net = hypot(ps(find(win,1,'last'),1) - ps(find(win,1),1), ps(find(win,1,'last'),2) - ps(find(win,1),2));
    pm.tortuosity = sum(seg) / max(net, eps);
    j = find(-vel(:,3) >= 0.9*pm.descentSpeed, 1);
    if ~isempty(j); pm.transientDrop = zDrop(j); end

    if pm.glideDeg >= opts.glideDeg;     pm.pathClass = "glide";
    elseif h.valid;                      pm.pathClass = "helix";
    elseif h.rms < opts.straightRms;     pm.pathClass = "straight";
    else;                                pm.pathClass = "irregular";
    end
    pm.valid = true;
end
