%% Flat-equivalence regression for the 3D-shape model
% Locks the invariant that the non-planar model (setupSeedShape3D + seed6DOFODE3D)
% reduces EXACTLY to the planar model (setupSeedShapeAndMass + seed6DOFODE) for a
% FLAT seed. Run this after any change to the 3D physics: if it still passes, the
% change did not disturb flat-plate dynamics.
%
% SCOPE NOTE (Sep-2026 audit). The two models are no longer identical out of the
% box, by design: the 3D RHS dropped the invented whole-seed terms (Tx, Ty, the
% span torque with its migrating CoP and reduced-frequency attenuation) while the
% planar RHS keeps them as the frozen reference. So equivalence is now asserted
% over the SHARED CORE -- Andersen-Pesavento-Wang sectional coefficients + strip
% theory + rigid-body dynamics + added mass -- by configuring the planar model
% down to it. That is a stronger statement than it looks: it says the cuts removed
% exactly the intended terms and touched nothing else.
%
% Checks:
%   1. The 'shape3d' flat seed has identical strips/mass to the planar seed, plus
%      y=0 and identity frames, and satisfies the seed-model contract. It must
%      also NOT carry the switches the 3D RHS no longer honours.
%   2. Core-only: seed6DOFODE3D(shape3d flat) matches seed6DOFODE(planar) exactly
%      when the planar extras are switched off.
%   3. The same for a planar seedParams fed through the 3D RHS (identity fallback).
%   4. The same for a moving-nut (time-varying mass) seed.
%   5. Anti-regression: with the planar extras switched back ON the two models MUST
%      now differ. If this difference ever returns to zero, a cut term has crept
%      back into the 3D RHS.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

tol = 1e-12;

% --- Base seed -----------------------------------------------------------
cfg.spanLength=0.050; cfg.chordLength=0.015; cfg.thickness=0.002; cfg.bulkDensity=65;
cfg.numStrips=10; cfg.tSamples=0; cfg.nutMass=75e-6; cfg.rhoFluid=1.225; cfg.g=9.81;
xh=cfg.spanLength/2; yh=cfg.chordLength/2;
baseBsp.seedShape=polyshape([-xh,xh,xh,-xh],[-yh,-yh,yh,yh]);
baseBsp.seedDensity=cfg.bulkDensity*cfg.thickness; baseBsp.seedThickness=cfg.thickness;
baseBsp.numStrips=cfg.numStrips;

% --- Configs -------------------------------------------------------------
% CORE: planar stripped down to what the 3D model still has. The span force is
% switched off on BOTH sides (it survives in the 3D model, but its moment is now
% computed differently, so it is not part of the shared core).
cfgPcore = cfg;
cfgPcore.enableSpanForce = false;  cfgPcore.enableTxDamping = false;
cfgPcore.enableNormalSpinDamping = false;
cfg3core = cfg;  cfg3core.shapeModel = 'shape3d';  cfg3core.enableSpanForce = false;

% FULL: each model at its own defaults (planar keeps the invented terms).
cfgPfull = cfg;
cfg3full = cfg;  cfg3full.shapeModel = 'shape3d';

%% 1. Build both, compare geometry/mass, validate contract
bsp = baseBsp; bsp.tSamples=0; bsp.nutPos_t=[0.4*cfg.chordLength;0;0.3*cfg.spanLength]; bsp.nutMass_t=cfg.nutMass;
spP = buildSeedParams(bsp, cfgPcore);      % planar, core-only config
sp3 = buildSeedParams(bsp, cfg3core);      % shape3d (flat), core-only config

iP = validateSeedParams(spP);   i3 = validateSeedParams(sp3);
fprintf('contract: planar model=%s, shape3d model=%s\n', iP.model, i3.model);

geomOK = isequal(spP.strips.chord, sp3.strips.chord) && isequal(spP.strips.dz, sp3.strips.dz) && ...
         isequal(spP.strips.xgc_body, sp3.strips.xgc_body) && isequal(spP.strips.zgc_body, sp3.strips.zgc_body) && ...
         isequal(spP.strips.liftMult, sp3.strips.liftMult) && isequal(spP.strips.dragMult, sp3.strips.dragMult) && ...
         isequal(spP.massParams.com_t, sp3.massParams.com_t) && isequal(spP.massParams.I_G_t, sp3.massParams.I_G_t) && ...
         isequal(spP.massParams.I_G_dot_t, sp3.massParams.I_G_dot_t) && isequal(spP.massParams.M_total, sp3.massParams.M_total);
framesFlat = all(sp3.strips.ygc_body==0) && isequal(sp3.strips.chordDir, repmat([1;0;0],1,i3.numStrips)) && ...
             isequal(sp3.strips.normalDir, repmat([0;1;0],1,i3.numStrips)) && isequal(sp3.strips.spanDir, repmat([0;0;1],1,i3.numStrips));

% The 3D builder must strip the switches whose terms were cut, so a stale config
% is visible rather than silently inert.
deadSwitches = {'enableSpanTorque','enableSpanCOPMigration', ...
                'enableSpanTorqueAttenuation','enableTxDamping','enableNormalSpinDamping'};
stripped = ~any(isfield(sp3, deadSwitches));
fprintf('geometry/mass identical to planar: %d   |   flat frames (y=0, identity): %d\n', geomOK, framesFlat);
fprintf('shape3d struct free of the cut switches: %d\n', stripped);

%% 2-3. RHS match on random states (fixed nut, CORE config)
rng(7); maxA=0; maxB=0;
for k=1:400
    x=[(rand(3,1)-0.5)*2; randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100]; x(4:7)=x(4:7)/norm(x(4:7));
    maxA = max(maxA, max(abs(seed6DOFODE3D(0,x,sp3) - seed6DOFODE(0,x,spP))));   % shape3d flat vs planar
    maxB = max(maxB, max(abs(seed6DOFODE3D(0,x,spP) - seed6DOFODE(0,x,spP))));   % planar seed thru 3D RHS
end
fprintf('\nCORE match (fixed nut, 400 states):\n');
fprintf('  seed6DOFODE3D(shape3d) vs seed6DOFODE(planar) : %.3e\n', maxA);
fprintf('  seed6DOFODE3D(planar)  vs seed6DOFODE(planar) : %.3e\n', maxB);

%% 4. RHS match on a moving-nut (time-varying mass) seed, CORE config
tS=0:0.02:2; bspM=bsp; bspM.tSamples=tS;
bspM.nutPos_t=[linspace(0,cfg.chordLength,numel(tS)); zeros(1,numel(tS)); linspace(0,cfg.spanLength,numel(tS))];
bspM.nutMass_t=cfg.nutMass*ones(size(tS));
spPm=buildSeedParams(bspM,cfgPcore); sp3m=buildSeedParams(bspM,cfg3core);
maxC=0;
for k=1:200
    x=[(rand(3,1)-0.5)*2; randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100]; x(4:7)=x(4:7)/norm(x(4:7)); t=rand*2;
    maxC=max(maxC, max(abs(seed6DOFODE3D(t,x,sp3m) - seed6DOFODE(t,x,spPm))));
end
fprintf('CORE match (moving nut, 200 states): %.3e\n', maxC);

%% 5. Anti-regression: the cut terms must actually be gone
% At full defaults the planar model applies Tx, Ty and the span torque; the 3D
% model applies none of them. The RHS difference is therefore the size of the
% terms that were cut -- it MUST be non-zero, or a cut has silently reverted.
spPf = buildSeedParams(bsp, cfgPfull);  sp3f = buildSeedParams(bsp, cfg3full);
maxD = 0;
for k=1:400
    x=[(rand(3,1)-0.5)*2; randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100]; x(4:7)=x(4:7)/norm(x(4:7));
    maxD = max(maxD, max(abs(seed6DOFODE3D(0,x,sp3f) - seed6DOFODE(0,x,spPf))));
end
cutsLive = maxD > 1e-6;
fprintf('\nFULL-default difference (= the cut terms): %.3e   -> cuts still in place: %d\n', maxD, cutsLive);

%% Verdict
pass = geomOK && framesFlat && stripped && iP.ok && i3.ok && ...
       maxA<tol && maxB<tol && maxC<tol && cutsLive;
if pass
    fprintf('\nPASS: the 3D model matches the planar CORE bit-for-bit (< %.0e), carries none\n', tol);
    fprintf('      of the cut switches, and still differs from the full planar model.\n');
else
    fprintf(2,'\nFAIL: flat-equivalence broken (see values above).\n');
end
