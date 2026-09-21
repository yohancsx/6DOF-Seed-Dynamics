%% Regression: curvature-correct mass properties (setupSeedShape3D)
% Locks the 3D mass/inertia derivation for a seed whose strips leave the body
% x-z plane. Run after any change to the shape builder's mass code.
%
% Checks:
%   A. ANALYTIC -- a uniformly-sliced flat rectangular plate must reproduce the
%      closed-form thin-plate tensor (1/12)m*diag(b^2, a^2+b^2, a^2). This
%      exercises the local-inertia formula, the perpendicular-axis theorem, and
%      the per-strip parallel-axis summation (which is exact for a uniform slice).
%   B. REDUCTION -- routing a FLAT seed through the curved-mass path must match
%      the planar builder's CoM / inertia / total mass.
%   C. CURVED SANITY -- a symmetric bowl: CoM stays on the symmetry axes but
%      shifts in +y (matching an independent hand calculation), I_G stays
%      symmetric positive-definite, and total mass is conserved.
%   D. TOGGLES -- each model stamps its own switch set (planar: span force on,
%      added-mass rate off; shape3d: edge drag + added-mass rate on, no span-force
%      switch), a cfg override wins, and the general per-strip added-mass form
%      reduces exactly to the flat form on a flat seed.
%
% A "flat" curvature profile (@(s) 0*s) is used for A and B: curvature is
% PRESENT (so the curved-mass path runs) but zero (so the geometry is flat).

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

tolRel = 1e-10;   tolAbs = 1e-12;

% --- Base seed -----------------------------------------------------------
cfg.spanLength=0.050; cfg.chordLength=0.015; cfg.thickness=0.002; cfg.bulkDensity=65;
cfg.numStrips=12; cfg.tSamples=0; cfg.nutMass=75e-6; cfg.rhoFluid=1.225; cfg.g=9.81;
a = cfg.chordLength;  b = cfg.spanLength;  xh = b/2;  yh = a/2;  hs = b/2;
rhoA = cfg.bulkDensity * cfg.thickness;
base.seedShape   = polyshape([-xh,xh,xh,-xh], [-yh,-yh,yh,yh]);
base.seedDensity = rhoA;  base.seedThickness = cfg.thickness;  base.numStrips = cfg.numStrips;
base.tSamples = 0;  base.nutPos_t = [0;0;0];  base.nutMass_t = cfg.nutMass;
cfg3 = cfg;  cfg3.shapeModel = 'shape3d';
flatCurv = @(s) 0*s;

%% A. Analytic: flat plate (no nut) vs the closed-form thin-plate tensor
bA = base;  bA.nutMass_t = 0;  bA.curvature = flatCurv;
spA = buildSeedParams(bA, cfg3);
m         = rhoA * a * b;
Ianalytic = (1/12) * m * diag([b^2, a^2 + b^2, a^2]);
IA        = spA.massParams.I_G_t(:,:,1);
errA      = max(abs(IA - Ianalytic), [], 'all') / max(abs(Ianalytic), [], 'all');
offA      = max(abs(IA - diag(diag(IA))), [], 'all');
okA = errA < tolRel && offA < tolAbs && abs(spA.massParams.M_total - m) < tolAbs;
fprintf('A) flat plate vs closed form : rel err %.2e, max offdiag %.2e, mass err %.2e  -> %s\n', ...
        errA, offA, abs(spA.massParams.M_total-m), pf(okA));

%% B. Reduction: curved-mass path on a flat seed vs the planar builder
spP = buildSeedParams(base, cfg);
bB  = base;  bB.curvature = flatCurv;   spB = buildSeedParams(bB, cfg3);
dCom = max(abs(spP.massParams.com_t  - spB.massParams.com_t),  [], 'all');
dI   = max(abs(spP.massParams.I_G_t  - spB.massParams.I_G_t),  [], 'all');
relI = dI / max(abs(spP.massParams.I_G_t), [], 'all');
dM   = abs(spP.massParams.M_total - spB.massParams.M_total);
okB  = dCom < tolAbs && relI < tolRel && dM < tolAbs;
fprintf('B) curved path on flat seed  : dCoM %.2e m, dI rel %.2e, dM %.2e            -> %s\n', ...
        dCom, relI, dM, pf(okB));

%% C. Curved sanity: symmetric bowl
tipPhi = deg2rad(35);
bC = base;  bC.curvature = @(s) (tipPhi/hs) * s;
spC = buildSeedParams(bC, cfg3);
com = spC.massParams.com_t(:,1);   IC = spC.massParams.I_G_t(:,:,1);
mi  = rhoA * spC.strips.area;      mw = sum(mi);
comY_hand = sum(mi .* spC.strips.ygc_body) / (mw + cfg.nutMass);   % nut at origin
symErr = max(abs(IC - IC.'), [], 'all');
okC = abs(com(1)) < tolAbs && abs(com(3)) < tolAbs && com(2) > 0 ...
      && abs(com(2) - comY_hand) < tolAbs && symErr < tolAbs && all(eig(IC) > 0) ...
      && abs(spC.massParams.M_total - spB.massParams.M_total) < tolAbs;
fprintf('C) bowl: CoM=[%.1e %.3e %.1e], hand-calc dy %.1e, symErr %.1e, eig>0 %d  -> %s\n', ...
        com(1), com(2), com(3), abs(com(2)-comY_hand), symErr, all(eig(IC)>0), pf(okC));
fprintf('   (curvature changes the inertia by %.1f%% vs flat)\n', ...
        100*max(abs(IC - spB.massParams.I_G_t(:,:,1)),[],'all')/max(abs(spB.massParams.I_G_t(:,:,1)),[],'all'));

%% D. Physics toggles: per-model switch sets + added-mass reduction
% Each model must stamp its OWN switch set. Planar keeps the span force and has
% the added-mass rate OFF (frozen reference). shape3d -- flat, twisted or curved
% alike -- has edge drag and the added-mass rate ON, and carries NO span-force
% switch at all (retired in phase 4; with it went the old curved-seed special
% case). A cfg override must win. The general added-mass form must still reduce
% EXACTLY to the flat form on a flat seed -- including with an out-of-plane nut,
% since the (d x n) term drops any offset parallel to the strip normal.
bT = base;  bT.twist = @(z) (deg2rad(20)/hs) * z;   spT = buildSeedParams(bT, cfg3);
cfgOv = cfg3;  cfgOv.enableEdgeDrag = false;
spCov = buildSeedParams(bC, cfgOv);
is3dSet = @(s) isfield(s,'enableEdgeDrag') && s.enableEdgeDrag ...
               && s.enableAddedMassRate && ~isfield(s,'enableSpanForce');
okD1 = spP.enableSpanForce && ~spP.enableAddedMassRate && ~isfield(spP,'enableEdgeDrag') ...
       && is3dSet(spT) && is3dSet(spC) && ~spCov.enableEdgeDrag;

spF = buildSeedParams(base, cfg3);                 % flat shape3d (identity frames)
bN  = base;  bN.nutPos_t = [0; -a; 0];   spN = buildSeedParams(bN, cfg3);  % out-of-plane nut
flatSeeds = {spF, spN};
dAM = 0;
for q = 1 : numel(flatSeeds)
    sq = flatSeeds{q};   cm = sq.massParams.com_t(:,1);
    s0 = sq;  s0.enableAddedMass3D = false;  [At0, Ar0] = getAddedMass(s0, cm);
    s1 = sq;  s1.enableAddedMass3D = true;   [At1, Ar1] = getAddedMass(s1, cm);
    dAM = max([dAM, max(abs(At1-At0),[],'all')/max(abs(At0),[],'all'), ...
                     max(abs(Ar1-Ar0),[],'all')/max(abs(Ar0),[],'all')]);
end
okD2 = dAM < tolRel;
okD = okD1 && okD2;
fprintf(['D) toggles: planar [spanForce %d, AdotRate %d]  shape3d twist/curved ' ...
         '[edgeDrag %d/%d, AdotRate %d/%d]  override edgeDrag->%d, addedMass reduction %.1e -> %s\n'], ...
        spP.enableSpanForce, spP.enableAddedMassRate, spT.enableEdgeDrag, spC.enableEdgeDrag, ...
        spT.enableAddedMassRate, spC.enableAddedMassRate, spCov.enableEdgeDrag, dAM, pf(okD));

%% Verdict
if okA && okB && okC && okD
    fprintf('\nPASS: curved mass/inertia is analytic-exact, reduces to the planar builder,\n');
    fprintf('      and behaves correctly for a curved seed.\n');
else
    fprintf(2, '\nFAIL: see the checks above.\n');
end

function s = pf(ok);  if ok; s = 'PASS'; else; s = 'FAIL'; end;  end
