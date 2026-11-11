function res = testTrackedSeedLoader()
% TESTTRACKEDSEEDLOADER  The tracked-drop import (loadTrackedSeeds) and the shared
%   position-only metrics (computePathMetrics) used to compare simulation with the
%   real drops (runRealSeedSet3D, suite stage G).
%
% CASES
%   T1  Import: t = 0 at release and uniform 1/60 s steps; zDrop starts at 0; z is
%       up (every drop falls); the partial flag marks exactly the tracks that end
%       early (seeds that glided out of view).
%   T2  Synthetic helix + drift + tracking noise at the level of the real tracks:
%       descent speed, drift, helix radius, period and sense are recovered.
%   T3  Synthetic glide (50 deg) with noise: glide angle recovered, class 'glide'.
%   T4  Synthetic straight descent with noise: no helix is claimed.
%   T5  Camera sampling: a simulated real-seed drop measured from its dense output
%       and from 60 fps samples gives the same metrics.
% Skipped with a WARN if the tracked-data file is absent.
%
% OUTPUT
%   res : struct array (.name, .pass, .detail).

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  Pth = seedPaths(here);
    res = struct('name', {}, 'pass', {}, 'detail', {});
    if ~isfile(Pth.seedTracks)
        res = item('tracked-seed loader', true, ['WARN skipped: ' Pth.seedTracks ' not found']);
        report(res);  return
    end

    % ============================== T1 ===========================================
    tr = loadTrackedSeeds(Pth.seedTracks);
    ok = true;  maxStep = 0;  minFall = Inf;
    for k = 1:numel(tr)
        d = diff(tr(k).t) * 60;
        ok = ok && tr(k).t(1) == 0 && tr(k).zDrop(1) == 0 && all(d >= 1 - 1e-9);
        maxStep = max(maxStep, max(abs(d(d < 1.5) - 1)));
        minFall = min(minFall, tr(k).fallen);
    end
    partial = sort([tr([tr.partial]).seedNum]);
    noise = trackNoise(tr);
    res(end+1) = item('T1 import: time base, z up, partial tracks', ok && maxStep < 1e-9 && minFall > 3, ...
        sprintf('%d drops; fallen >= %.1f m; partial (left the view): %s; position noise ~%.3f m', ...
                numel(tr), minFall, mat2str(partial), noise));

    % ============================== T2 / T3 / T4 =================================
    rng(41);  fps = 60;  t = (0:1/fps:4).';
    Vd = 2.1;  R = 0.25;  Tp = 1.2;  drift = 0.2;  w = 2*pi/Tp;
    pos = [drift*t + R*cos(w*t), R*sin(w*t), 10 - Vd*t];
    pos = pos + noise*randn(size(pos));
    p = computePathMetrics(t, pos, struct('window', [1.5, 8]));
    h = p.helix;
    e2 = [abs(p.descentSpeed/Vd - 1), abs(p.driftSpeed/drift - 1), abs(h.radius/R - 1), abs(h.period/Tp - 1)];
    res(end+1) = item('T2 synthetic helix recovered', all(e2 < [0.02 0.15 0.15 0.10]) && h.turns > 0 && ...
        p.pathClass == "helix", sprintf('rel err Vd %.3f, drift %.3f, R %.3f, period %.3f; turns %+.2f; %s', ...
        e2, h.turns, p.pathClass));

    gl = 50;  vh = Vd*tand(gl);
    pos = [vh*t, zeros(size(t)), 10 - Vd*t] + noise*randn(numel(t), 3);
    p = computePathMetrics(t, pos, struct('window', [1.5, 8]));
    res(end+1) = item('T3 synthetic glide recovered', abs(p.glideDeg - gl) < 2 && p.pathClass == "glide", ...
        sprintf('glide %.1f deg (true %d), class %s', p.glideDeg, gl, p.pathClass));

    pos = [zeros(numel(t), 2), 10 - Vd*t] + noise*randn(numel(t), 3);
    p = computePathMetrics(t, pos, struct('window', [1.5, 8]));
    res(end+1) = item('T4 straight descent: no helix claimed', ~p.helix.valid && p.pathClass ~= "helix", ...
        sprintf('class %s, helix rms %.3f m, explained %.2f', p.pathClass, p.helix.rms, p.helix.explained));

    % ============================== T5 ===========================================
    try
        s  = loadRealSeeds(Pth.seedMassModel, struct('seeds', 10));
        sp = buildSeedParams(s.bsp, struct('rhoFluid', 1.225, 'g', 9.81, 'enableAddedMass3D', true));
        x0 = [0; 6; 0; axisAngleToQuat([0;0;1], deg2rad(135)); zeros(6,1)];
        opt = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Events', @ground);
        [ts, xs] = ode45(@(tt,xx) seed6DOFODE3D(tt,xx,sp), [0 20], x0, opt);
        toLab = @(X) [X(:,1), -X(:,3), X(:,2)];
        tD = (0:1/600:ts(end)).';   t6 = (0:1/60:ts(end)).';
        pD = computePathMetrics(tD, toLab(interp1(ts, xs(:,1:3), tD)), struct('window', [1.5, 5.5]));
        p6 = computePathMetrics(t6, toLab(interp1(ts, xs(:,1:3), t6)), struct('window', [1.5, 5.5]));
        e5 = [abs(p6.descentSpeed/pD.descentSpeed - 1), abs(p6.glideDeg - pD.glideDeg), ...
              abs(p6.helix.rms/pD.helix.rms - 1)];
        res(end+1) = item('T5 60 fps sampling = dense output', all(e5 < [0.01 1 0.05]), sprintf( ...
            'Vd rel %.1e, glide %.2f deg, helix rms rel %.3f (seed 10, 135 deg release)', e5));
    catch err
        res(end+1) = item('T5 60 fps sampling = dense output', false, err.message);
    end

    report(res);
end

function [value, isterminal, direction] = ground(~, x)
    value = x(2);  isterminal = 1;  direction = -1;
end
function s = trackNoise(tr)
% Position noise of the tracking, from the scatter of third differences of z
% (smooth motion contributes almost nothing at 60 fps): var(d3) = 20 sigma^2.
    v = [];
    for k = 1:numel(tr); v = [v; diff(tr(k).pos(:,3), 3)]; end %#ok<AGROW>
    s = std(v) / sqrt(20);
end
function report(res)
    for k = 1:numel(res)
        if ~res(k).pass; tag = 'FAIL'; elseif startsWith(res(k).detail, 'WARN'); tag = 'WARN'; else; tag = 'PASS'; end
        fprintf('  %-44s %s  %s\n', res(k).name, tag, res(k).detail);
    end
end
function r = item(name, pass, detail)
    r = struct('name', name, 'pass', logical(pass), 'detail', char(detail));
end
