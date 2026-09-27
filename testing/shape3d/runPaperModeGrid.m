%% Phase map over the Hou et al. (2025) parameter window -- our model vs their Fig. 2a
% A SEPARATE grid from testing/planar/runSeedModeGrid.m, which is untouched so that
% every earlier run stays comparable. Differences, all deliberate:
%   * their plate and mass ratios (paperSeedConfig), not our 15 mm seed;
%   * their axes -- CoM position as (x_c/a, y_c/b), spanwise and chordwise,
%     swept only over the window their weights can actually reach
%     (x_c/a <= 0.416, y_c/b <= 0.208);
%   * their release (long axis horizontal, face vertical, from rest);
%   * flat shape3d via seedRHS -- they excluded camber and twist;
%   * BOTH classifiers reported per cell, theirs and ours, plus a cross-tab.
%
% Scored against paperModeReference, a hand-digitised reading of their figure --
% see that file for how rough that is. The point of this grid is the TOPOLOGY:
% which mode sits where, and whether the AR/FA and ST/CH boundaries land in about
% the right place.
%
% Figures and the .mat go OUTSIDE the repo (generated binaries stay out of git).

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));
outDir = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Paper Mode Grid";
if ~exist(outDir,'dir'); mkdir(outDir); end

%% 0. Configuration -- EDIT HERE
[cfg, bsp, toNutPos, info] = paperSeedConfig();

% PHYSICS: left at the shape3d builder's defaults. The Tier 1 screen picked R+M
% (the two halves of Kirchhoff's pair, both in APW eqs. 6.1-6.3) and they are now
% the defaults, so there is nothing to override here -- pinning them again would
% just let this script drift from the model. LEV is off by default: it was
% near-neutral with M on and actively harmful without it.
% Read back whatever the builder actually stamped, so the run records the physics
% it was run with rather than a hand-copied guess at it.
bspSample = bsp;  bspSample.nutPos_t = [0;0;0];
spSample  = buildSeedParams(bspSample, cfg);
cfg.enableAddedMassRate   = spSample.enableAddedMassRate;
cfg.enableAddedMassMoment = spSample.enableAddedMassMoment;
cfg.enableLEV             = spSample.enableLEV;

nX = 13;   nY = 8;                   % grid resolution (13 x 8 = 104 runs)
xaGrid = linspace(0.01, 0.41, nX);   % x_c/a : spanwise  CoM offset
ybGrid = linspace(0.01, 0.20, nY);   % y_c/b : chordwise CoM offset

mopts = cfg.metricOpts;
mopts.refLength = info.L;   mopts.sigma = info.sigma;
mopts.g         = cfg.g;    mopts.rhoFluid = cfg.rhoFluid;

%% 1. Run the grid
paperMode = strings(nY, nX);   ourMode = strings(nY, nX);
refMode   = strings(nY, nX);   Vd      = nan(nY, nX);
theta     = nan(nY, nX);       thetaSd = nan(nY, nX);
fprintf('Paper mode grid: %d x %d = %d runs (R%d M%d L%d)\n', nX, nY, nX*nY, ...
        cfg.enableAddedMassRate, cfg.enableAddedMassMoment, cfg.enableLEV);
tStart = tic;
for iy = 1:nY
    for ix = 1:nX
        refMode(iy,ix) = string(paperModeReference(xaGrid(ix), ybGrid(iy)));
        try
            [~, r] = evalc(['runSingleMode('''', toNutPos(xaGrid(ix), ybGrid(iy)), ' ...
                            'cfg.releaseQuat, cfg.releaseOmega, cfg, bsp)']);
            if any(~isfinite(r.x(:)))
                paperMode(iy,ix) = "failed";  ourMode(iy,ix) = "failed";
            else
                pm = computePaperMetrics(r.t, r.x, mopts);
                paperMode(iy,ix) = string(classifyPaperMode(pm));
                ourMode(iy,ix)   = string(r.mode);
                Vd(iy,ix)        = pm.descentSpeed;
                theta(iy,ix)     = pm.spanAxisTiltDeg;
                thetaSd(iy,ix)   = pm.spanAxisTiltStd;
            end
        catch ME
            paperMode(iy,ix) = "failed";   ourMode(iy,ix) = "failed";
            fprintf('  cell (%.2f, %.2f) failed: %s\n', xaGrid(ix), ybGrid(iy), ME.message);
        end
    end
    fprintf('  row %d/%d (y_c/b = %.3f) done  [%.0f s]\n', iy, nY, ybGrid(iy), toc(tStart));
end

%% 2. SAVE FIRST -- before any reporting
% The grid is the expensive part (~25 min). Reporting is cheap and is where the
% bugs are, so the results go to disk immediately: a reporting error must never
% cost a re-run. reportPaperModeGrid.m reads this file and can be re-run alone.
gridFile = fullfile(outDir, 'paper_mode_grid.mat');
save(gridFile, 'xaGrid','ybGrid','paperMode','ourMode','refMode','Vd', ...
     'theta','thetaSd','cfg','info','outDir');
fprintf('\nGrid finished in %.0f s and saved to %s\n', toc(tStart), gridFile);

%% 3. Report (separate script, so it can be iterated without re-running)
reportPaperModeGrid(gridFile);
