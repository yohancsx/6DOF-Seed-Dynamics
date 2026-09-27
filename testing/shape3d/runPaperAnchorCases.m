%% Anchor cases: the five configurations Hou et al. (2025) publish kinematics for
% Drops OUR model at the exact (x_c/a, y_c/b) of each representative case in Hou,
% Zhang, Li, Jia & Huang, Commun. Eng. 4, 129 (2025), on THEIR plate (see
% paperSeedConfig) and at THEIR release condition, then reports our kinematics
% next to their published numbers.
%
% This is the phase-5 measuring stick. Two uses:
%   1. NOW: quantify the baseline mismatch with the current shape3d defaults
%      (edge drag ON, added-mass rate OFF, LEV OFF) -- no new physics.
%   2. NEXT: score Tier 1 switch combinations cheaply (5 runs each) before
%      committing to a full grid.
% It also supplies the numbers that the paper-mode classifier's thresholds will be
% read off, which is why it comes BEFORE that classifier exists.
%
% THEIR PUBLISHED TARGETS (Results > "Kinematics of typical flight patterns"):
%   AR  (0.39, 0.17)  theta = -11.4 +/- 4.2 deg, no tumbling, LOWEST descent speed
%   CST (0.06, 0.04)  theta = -38.2 +/- 2.3 deg, ~7 tumbles per revolution, R ~ 1.25 L
%   SST (0.10, 0.10)  theta = -38.3 +/- 13.0 deg, R ~ 1.89 L, tumbling REVERSES at
%                     each turning point
%   CH  (0.18, 0.08)  irregular, no sustained periodicity, large scatter
%   FA  (0.35, 0.08)  HIGHEST descent speed, large tilt, heavy side down
% Ordering they measure: V_d(AR) < V_d(CST) < V_d(SST) < V_d(FA), with CH scattered.
% Their scaling law (eq. 2) puts V_d / sqrt(sigma*g/rho) at O(1) for every mode.
%
% Sign note: their theta is negative, meaning the tip on the weighted side sits
% below horizontal. computePaperMetrics uses the same convention (+z = nut side),
% so the signs are directly comparable.
%
% TWO CONFOUNDS CHECKED AND RULED OUT (Sep 2026), so they need not be re-run:
%   * RELEASE ATTITUDE. Their psi = pi/2 is the self-rotation about the long axis,
%     but the paper does not say where psi = 0 points, so face-down and
%     face-vertical are both readable. Running every anchor both ways gives the
%     SAME mode label in all five, theta within 5 deg and V_d within 0.8 m/s
%     (SST moves most, which is fair -- it is the multistable one). The choice
%     does not drive these results.
%   * STRIP COUNT. 12 vs 24 strips: theta within 1 deg and V_d within 0.07 m/s on
%     every case. Converged.
%   What the strip check DID expose is a metric limitation: tumblesPerRev swung
%   200 -> 76 and 43 -> 96 on the erratic cases, because a wandering ground track
%   still admits a circle fit and the orbital rate is then meaningless. It is
%   trustworthy only where the cone is steady (CST 0.8 and FA 8.1 held across
%   every variation). Values from an unsteady cone are flagged '*' below.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

[cfg, bsp, toNutPos, info] = paperSeedConfig();

% --- OPTIONAL Tier 1 overrides -- leave commented for the baseline run -------
% cfg.enableAddedMassRate = true;
% cfg.enableLEV           = true;

mopts = cfg.metricOpts;
mopts.refLength = info.L;   mopts.sigma = info.sigma;
mopts.g         = cfg.g;    mopts.rhoFluid = cfg.rhoFluid;

%   name                 x_c/a  y_c/b  theta   thetaStd  tumb/rev  R/L
anchors = { 'AR  autorotation',   0.39,  0.17, -11.4,   4.2,      NaN,   NaN ; ...
            'CST continuous ST',  0.06,  0.04, -38.2,   2.3,      7.0,   1.25; ...
            'SST segmented ST',   0.10,  0.10, -38.3,  13.0,      NaN,   1.89; ...
            'CH  chaotic',        0.18,  0.08,   NaN,   NaN,      NaN,   NaN ; ...
            'FA  falling',        0.35,  0.08,   NaN,   NaN,      NaN,   NaN };

fprintf('\n=== Their plate, in our units (paperSeedConfig) ===\n');
fprintf('  span(L) %.0f mm x chord(W) %.0f mm, AR %.2f, t = %.2f mm, %d strips\n', ...
        info.L*1e3, info.W*1e3, info.L/info.W, info.thickness*1e3, cfg.numStrips);
fprintf('  plate %.1f mg + nut %.1f mg = %.1f mg,  mu = %.3f,  sigma = %.3f kg/m^2\n', ...
        info.plateMass*1e6, info.nutMass*1e6, info.totalMass*1e6, info.mu, info.sigma);
fprintf('  velocity scale sqrt(sigma g / rho) = %.3f m/s  (their V_d should be ~this)\n', info.Vscale);
fprintf('  release: long axis horizontal, face vertical (theta=0, psi=pi/2), from rest\n');
fprintf('  window caps: x_c/a <= %.3f, y_c/b <= %.3f\n\n', info.xcMax, info.ycMax);

%% Run
nA = size(anchors, 1);
res = cell(nA, 1);   pmAll = cell(nA, 1);   pm4All = cell(nA, 1);
for i = 1:nA
    nutPos = toNutPos(anchors{i,2}, anchors{i,3});
    r = runSingleMode(anchors{i,1}, nutPos, cfg.releaseQuat, cfg.releaseOmega, cfg, bsp);
    res{i}   = r;
    pmAll{i} = computePaperMetrics(r.t, r.x, mopts);

    % Same run read at THEIR observation window (a 4 m fall), to catch a mode that
    % only appears after the plate has fallen further than their rig allows.
    k4 = find(r.x(:,2) <= -cfg.dropHeight, 1);
    if ~isempty(k4) && k4 > 20
        pm4All{i} = computePaperMetrics(r.t(1:k4), r.x(1:k4,:), mopts);
    else
        pm4All{i} = [];
    end
end

%% Table 1 -- attitude and rotation vs their published kinematics
fprintf('\n=== Attitude and rotation (steady-state window) ===\n');
fprintf('%-19s %-5s %-5s | %-17s %-17s | %-9s %-7s | %-6s %-4s\n', ...
        'case','x_c/a','y_c/b','theta model (deg)','theta paper (deg)', ...
        'tumb/rev','paper','revers','ours');
fprintf('%s\n', repmat('-', 1, 112));
for i = 1:nA
    pm = pmAll{i};
    tgt = fmtPM(anchors{i,4}, anchors{i,5});
    tpr = fmtNum(anchors{i,6}, '%.1f');
    fprintf('%-19s %-5.2f %-5.2f | %-17s %-17s | %-9s %-7s | %-6d %-4s\n', ...
            anchors{i,1}, anchors{i,2}, anchors{i,3}, ...
            sprintf('%+.1f +/- %.1f', pm.spanAxisTiltDeg, pm.spanAxisTiltStd), tgt, ...
            [fmtNum(pm.tumblesPerRev, '%.2f') unreliable(pm)], tpr, ...
            pm.tumbleReversals, res{i}.mode);
end

%% Table 2 -- speed and track
fprintf('\n=== Descent speed and ground track ===\n');
fprintf('%-19s | %-7s %-10s | %-8s %-7s | %-9s %-6s | %-8s\n', ...
        'case','V_d','V_d/scale','R/L','paper','tumble w_z','conv','V_d @4m');
fprintf('%s\n', repmat('-', 1, 100));
for i = 1:nA
    pm = pmAll{i};   pm4 = pm4All{i};
    if isempty(pm4); v4 = NaN; else; v4 = pm4.descentSpeed; end
    fprintf('%-19s | %-7.2f %-10.2f | %-8s %-7s | %-9.2f %-6d | %-8s\n', ...
            anchors{i,1}, pm.descentSpeed, pm.VdScaled, ...
            fmtNum(pm.helixRadiusPerL, '%.2f'), fmtNum(anchors{i,7}, '%.2f'), ...
            pm.tumbleRate, pm.converged, fmtNum(v4, '%.2f'));
end

%% Checks the paper states outright
fprintf('\n=== Ordering and scale checks ===\n');
Vd = cellfun(@(p) p.descentSpeed, pmAll);
[~, ordIdx] = sort(Vd);
fprintf('  descent-speed order (slowest first): %s\n', ...
        strjoin(cellfun(@(k) strtrim(anchors{k,1}(1:4)), num2cell(ordIdx).', 'uni', 0), ' < '));
fprintf('  they report:                         AR < CST < SST < FA (CH scattered)\n');
fprintf('  V_d / sqrt(sigma g / rho): min %.2f, max %.2f  (their prefactor is O(1))\n', ...
        min(cellfun(@(p) p.VdScaled, pmAll)), max(cellfun(@(p) p.VdScaled, pmAll)));
okAR = Vd(1) == min(Vd);
okFA = Vd(5) == max(Vd);
fprintf('  AR slowest: %d    FA fastest: %d\n', okAR, okFA);

outFile = fullfile(tempdir, 'paperAnchorCases.mat');
save(outFile, 'anchors', 'pmAll', 'pm4All', 'cfg', 'info');
fprintf('\nMetrics saved to %s\n', outFile);

% =========================================================================
function s = fmtPM(v, sd)
    if isnan(v);  s = '(not reported)';  else;  s = sprintf('%+.1f +/- %.1f', v, sd);  end
end
function s = fmtNum(v, f)
    if isnan(v) || isinf(v);  s = '--';  else;  s = sprintf(f, v);  end
end
function s = unreliable(pm)
% '*' marks an orbital-rate-derived number taken from a cone that is not steady,
% where the circle fit -- and so tumblesPerRev -- carries no real information.
    if pm.spanAxisTiltStd > 10;  s = '*';  else;  s = '';  end
end
