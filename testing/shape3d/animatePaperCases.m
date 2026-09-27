%% Animate representative paper-window cases, to check the labels BY EYE
% The paper-mode classifier decides every conclusion in the phase-5 comparison,
% and testPaperMetrics proves it is right on synthetic motion where the answer is
% known. This script closes the other half: it renders the ACTUAL simulated
% trajectories so the label can be checked against what the seed visibly does.
%
% For each case it prints the metrics the classifier used, then writes a video.
% What to look for, per label:
%   AR   revolves about a vertical axis, long axis at a steady shallow tilt, and
%        does NOT flip over. (Watch for the trap: a tilted revolution has a
%        constant nonzero omega_z without tumbling at all -- see testPaperMetrics
%        case A. If it reads AR but you SEE it flipping end-over-end, that is the
%        bug this animation exists to catch.)
%   CST  flips end-over-end about its LONG axis, continuously, one direction,
%        while the whole path spirals down a vertical axis at a steady cone.
%   SST  the same, but the flip direction REVERSES at turning points.
%   FA   falls at a large, steady tilt without sustained rotation of either kind.
%   CH   no sustained pattern; attitude swings.
%
% Videos go OUTSIDE the repo (generated binaries stay out of git).

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));
outDir = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Paper Case Animations';
if ~exist(outDir,'dir'); mkdir(outDir); end

[cfg, bsp, toNutPos, info] = paperSeedConfig();
mopts = cfg.metricOpts;
mopts.refLength = info.L;   mopts.sigma = info.sigma;
mopts.g         = cfg.g;    mopts.rhoFluid = cfg.rhoFluid;

% The five published anchors, plus two grid cells that sit on either side of the
% bifurcation our map produces (see the phase-5 notes: the map is effectively 1D
% in x_c/a, with a flat state to the left and a coned autorotation to the right).
cases = { 'anchor_AR',        0.39, 0.17 ; ...
          'anchor_CST',       0.06, 0.04 ; ...
          'anchor_SST',       0.10, 0.10 ; ...
          'anchor_CH',        0.18, 0.08 ; ...
          'anchor_FA',        0.35, 0.08 ; ...
          'grid_left_of_jump',0.14, 0.064; ...
          'grid_right_of_jump',0.24,0.064 };

fprintf('\n%-20s %-5s %-5s | %-5s | %-9s %-8s %-8s %-7s %-6s\n', ...
        'case','x/a','y/b','label','theta','psi''','phi''','turns','revs');
fprintf('%s\n', repmat('-', 1, 92));

for i = 1:size(cases,1)
    nutPos = toNutPos(cases{i,2}, cases{i,3});
    [~, r] = evalc('runSingleMode('''', nutPos, cfg.releaseQuat, cfg.releaseOmega, cfg, bsp)');
    pm     = computePaperMetrics(r.t, r.x, mopts);
    [lab, linfo] = classifyPaperMode(pm);

    fprintf('%-20s %-5.2f %-5.2f | %-5s | %+7.1f%s %-8.1f %-8.2f %-7.1f %-6d\n', ...
            cases{i,1}, cases{i,2}, cases{i,3}, lab, ...
            pm.spanAxisTiltDeg, sprintf('+/-%.1f', pm.spanAxisTiltStd), ...
            pm.selfRotationRate, pm.revolutionRate, pm.maxTurnsOneWay, pm.tumbleReversals);
    fprintf('    %s\n', linfo.reason);

    % --- animate -----------------------------------------------------------
    % Resample onto a uniform grid first: ode45's own steps cluster where the
    % dynamics are fast, which would make the playback speed misleading.
    tU = linspace(r.t(1), r.t(end), 600).';
    xU = interp1(r.t, r.x, tU);
    traj = struct('t', tU, 'comPos_world', xU(:,1:3).', 'quat', xU(:,4:7).', ...
                  'comVel_world', xU(:,8:10).', 'comOmega_body', xU(:,11:13).');

    opts = struct();
    opts.videoFile  = fullfile(outDir, sprintf('%s_%s.mp4', cases{i,1}, lab));
    opts.fps        = 25;
    opts.frameStep  = 2;
    opts.axisMode   = 'follow';
    opts.zoomFactor = 2.5;
    opts.showVels   = false;
    opts.figSize    = [720 540];
    opts.shapeOpts.showCom = true;
    opts.shapeOpts.showNut = true;
    animateSeed(r.seedParams, traj, opts);
end

fprintf('\nVideos written to %s\n', outDir);
fprintf(['Check each against the list at the top of this file. The one to watch is\n' ...
         'anchor_AR: it must revolve WITHOUT flipping end-over-end.\n']);
