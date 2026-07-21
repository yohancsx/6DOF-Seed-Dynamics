function m = computeTrajectoryMetrics(t, x, opts)
% COMPUTETRAJECTORYMETRICS  Reduce a seed trajectory to steady-state descriptors.
%
% Computes a small set of physically meaningful scalars that separate the
% biological flight modes (glide, dive, tumble, spiral, autorotation,
% parachute, flutter). All averages are taken over a STEADY-STATE WINDOW (the
% latter part of the run) so startup transients don't contaminate them.
%
% INPUTS
%   t    : Nx1 time vector (s).
%   x    : Nx13 state history [r(3) q(4) v(3) omega(3)] from seed6DOFODE (rows =
%          time). r/v are inertial (world), omega is body frame.
%   opts : (optional) struct:
%          .windowStartFrac - window is t >= this*t_end (default 0.5)
%          .convergeTol     - relative-change tolerance for "converged" (0.20)
%
% OUTPUT (struct m)
%   .descentSpeed     mean downward speed  (m/s, + = descending)
%   .horizontalSpeed  mean horizontal speed (m/s)
%   .glideRatio       net horizontal displacement / net height lost (dimensionless)
%   .straightness     net horizontal displacement / horizontal path length (0..1)
%   .helixRadius      radius of a circle fit to the horizontal ground track (m)
%   .helixValid       true if the circle fit is meaningful (curved path)
%   .spanwiseSpin     mean |omega_z| body  (rad/s)  -- end-over-end / tumbling
%   .normalSpin       mean |omega_y| body  (rad/s)  -- in-plane spin
%   .verticalSpinRate mean world-vertical spin (rad/s, signed -> handedness)
%   .verticalSpinMag  mean |world-vertical spin| (rad/s)
%   .coneAngleDeg     mean plate tilt from horizontal (deg; 0=broadside, 90=edge-on)
%   .tiltStd          std of that tilt over the window (deg) -- flutter proxy
%   .converged        true if descent speed & spin plateaued across the window
%   .windowStartFrac  the window fraction actually used

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'windowStartFrac'); opts.windowStartFrac = 0.5;  end
    if ~isfield(opts, 'convergeTol');     opts.convergeTol     = 0.20; end

    t = t(:);
    N = numel(t);
    r = x(:, 1:3);   q = x(:, 4:7);   v = x(:, 8:10);   w = x(:, 11:13);

    % --- Steady-state window (latter part of the run) ---------------------
    win = t >= opts.windowStartFrac * t(end);
    if nnz(win) < 10; win = true(N, 1); end   % too short -> use everything
    iw = find(win);
    nw = numel(iw);

    % --- Orientation-dependent quantities (need R per step) ---------------
    normalTilt = zeros(nw, 1);   % deg, plate normal from world vertical -> [0,90]
    vertSpin   = zeros(nw, 1);   % rad/s, angular velocity about world vertical (Y)
    for kk = 1:nw
        k = iw(kk);
        R = quatToRotm(q(k, :).');
        nWorld = R * [0; 1; 0];                          % plate normal in world
        normalTilt(kk) = acosd(min(1, abs(nWorld(2))));  % tilt from horizontal
        wWorld = R * w(k, :).';                          % body omega -> world
        vertSpin(kk) = wWorld(2);
    end

    % --- Translational descriptors (from the inertial velocity state) -----
    descentSpeed    = -mean(v(iw, 2));                          % + downward
    horizontalSpeed =  mean(hypot(v(iw, 1), v(iw, 3)));

    % --- Displacement-based glide ratio & path straightness ---------------
    X = r(iw, 1);  Y = r(iw, 2);  Z = r(iw, 3);
    horizDisp    = hypot(X(end) - X(1), Z(end) - Z(1));
    vertDrop     = Y(1) - Y(end);
    glideRatio   = horizDisp / max(vertDrop, eps);
    horizPath    = sum(hypot(diff(X), diff(Z)));
    straightness = horizDisp / max(horizPath, eps);

    % --- Helix radius (circle fit to the horizontal ground track) ---------
    [helixRadius, helixValid] = fitCircleRadius(X, Z);

    % --- Rotational descriptors -------------------------------------------
    spanwiseSpin     = mean(abs(w(iw, 3)));   % body z: end-over-end / tumbling
    normalSpin       = mean(abs(w(iw, 2)));   % body y: in-plane spin
    verticalSpinRate = mean(vertSpin);        % signed (handedness)
    verticalSpinMag  = mean(abs(vertSpin));

    % --- Attitude ---------------------------------------------------------
    coneAngleDeg = mean(normalTilt);          % 0 = flat/broadside, 90 = edge-on
    tiltStd      = std(normalTilt);           % oscillation amplitude (flutter)

    % --- Convergence: descent speed & spin steady across window halves ----
    half = floor(nw / 2);
    if half >= 3
        d1 = -mean(v(iw(1:half),     2));
        d2 = -mean(v(iw(half+1:end), 2));
        s1 =  mean(abs(vertSpin(1:half)));
        s2 =  mean(abs(vertSpin(half+1:end)));
        relD = abs(d2 - d1) / max(abs(d2), eps);
        relS = abs(s2 - s1) / max(abs(s2), 0.1);   % 0.1 rad/s floor vs noise blow-up
        converged = (relD < opts.convergeTol) && (relS < 3 * opts.convergeTol);
    else
        converged = false;
    end

    % --- Pack -------------------------------------------------------------
    m = struct('descentSpeed', descentSpeed, 'horizontalSpeed', horizontalSpeed, ...
               'glideRatio', glideRatio, 'straightness', straightness, ...
               'helixRadius', helixRadius, 'helixValid', helixValid, ...
               'spanwiseSpin', spanwiseSpin, 'normalSpin', normalSpin, ...
               'verticalSpinRate', verticalSpinRate, 'verticalSpinMag', verticalSpinMag, ...
               'coneAngleDeg', coneAngleDeg, 'tiltStd', tiltStd, ...
               'converged', converged, 'windowStartFrac', opts.windowStartFrac);
end


% =========================================================================
% LOCAL: algebraic (Kasa) circle fit; returns radius and a validity flag
% =========================================================================
function [R, valid] = fitCircleRadius(u, v)
    u = u(:);  v = v(:);
    if numel(u) < 5
        R = Inf; valid = false; return
    end
    A   = [u, v, ones(numel(u), 1)];
    b   = u.^2 + v.^2;
    sol = A \ b;
    uc  = sol(1) / 2;
    vc  = sol(2) / 2;
    R   = sqrt(max(sol(3) + uc^2 + vc^2, 0));
    % Valid only if the points actually lie near a circle of finite radius:
    % the spread of point-to-centre distances should be small relative to R.
    resid = std(hypot(u - uc, v - vc));
    valid = isfinite(R) && R > 0 && resid < 0.5 * R;
end
