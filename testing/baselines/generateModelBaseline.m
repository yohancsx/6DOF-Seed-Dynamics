%% Generate a model-baseline snapshot
% Records WHAT THE MODEL DOES at the current code/physics state, as a regression
% baseline for future changes (e.g. the planar -> full-3D work). Produces a
% timestamped, git-stamped folder under model_test_results/ containing:
%   - a manifest (git commit + dirty flag, MATLAB version, seed geometry, physics
%     switches, aero + classifier settings, and each stage's config/inputs),
%   - a COARSE flight-mode phase map (mode_grid/),
%   - the 3 moving-CoM scenarios (com_movement/),
%   - the sweep test suite (test_suite/).
% The git hash pins the exact code (incl. aero constants), so a snapshot is
% reproducible from source. Outputs are git-ignored (see model_test_results/).
%
% Run headless (matlab -batch) or interactively. Needs physics/, visualization/,
% and testing/helpers/ on the path. Uses parfor for the grid if available.

%% 0. Paths + snapshot folder
repoRoot = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(repoRoot,'physics'), fullfile(repoRoot,'physics','helpers'), ...
        fullfile(repoRoot,'physics','aero'), fullfile(repoRoot,'physics','mass'), ...
        fullfile(repoRoot,'physics3d'), fullfile(repoRoot,'visualization'), ...
        fullfile(repoRoot,'testing','helpers'));

githash = gitField(repoRoot, 'rev-parse --short HEAD');
gitdirty = ~isempty(strtrim(gitField(repoRoot, 'status --porcelain')));
stamp   = datestr(now, 'yyyy-mm-dd_HHMMSS');   %#ok<TNOW1,DATST>
snapDir = fullfile(repoRoot, 'model_test_results', sprintf('%s_%s', stamp, githash));
gridDir = fullfile(snapDir, 'mode_grid');
comDir  = fullfile(snapDir, 'com_movement');
tsDir   = fullfile(snapDir, 'test_suite');
cellfun(@(d) mkdir(d), {snapDir, gridDir, comDir, tsDir});
fprintf('Baseline snapshot: %s  (git %s%s)\n', snapDir, githash, ternary(gitdirty,'-dirty',''));

% --- Shared base seed (the working seed used everywhere) ------------------
base.spanLength = 0.050;  base.chordLength = 0.015;  base.thickness = 0.002;
base.bulkDensity = 65;    base.numStrips = 10;       base.nutMass = 75e-6;
base.rhoFluid = 1.225;    base.g = 9.81;
xh = base.spanLength/2;   yh = base.chordLength/2;
baseBsp.seedShape = polyshape([-xh,xh,xh,-xh],[-yh,-yh,yh,yh]);
baseBsp.seedDensity = base.bulkDensity*base.thickness;
baseBsp.seedThickness = base.thickness;  baseBsp.numStrips = base.numStrips;

%% 1. Manifest -- the exact physics/config so the runs can be replicated
% Read the switches off a freshly built seed (captures the code defaults).
bspSample = baseBsp;  bspSample.tSamples = 0;
bspSample.nutPos_t = [0;0;0];  bspSample.nutMass_t = base.nutMass;
spSample  = buildSeedParams(bspSample, struct('rhoFluid',base.rhoFluid,'g',base.g));
switches  = struct('enableSpanForce',spSample.enableSpanForce, ...
                   'enableSpanTorque',spSample.enableSpanTorque, ...
                   'enableSpanGeomVelocity',spSample.enableSpanGeomVelocity, ...
                   'enableSpanCOPMigration',spSample.enableSpanCOPMigration, ...
                   'enableSpanTorqueAttenuation',spSample.enableSpanTorqueAttenuation, ...
                   'enableTxDamping',spSample.enableTxDamping);
% Aero fingerprint at a few AoAs (default law; exact constants pinned by git hash).
aeroCoeffs   = computeAeroCoeffs(0, []);       % coeffs at alpha = 0
aeroCoeffs45 = computeAeroCoeffs(pi/4, []);    % coeffs at alpha = 45 deg
modeThresh   = defaultModeThresholds();

manifest = struct('stamp',stamp,'gitHash',githash,'gitDirty',gitdirty, ...
                  'matlab',version,'baseSeed',base,'switches',switches, ...
                  'aeroCoeffs',aeroCoeffs,'aeroCoeffs45',aeroCoeffs45,'modeThresholds',modeThresh);
save(fullfile(snapDir,'manifest.mat'),'-struct','manifest');
writeManifestTxt(fullfile(snapDir,'manifest.txt'), manifest);
fprintf('  wrote manifest\n');

%% 2. Coarse mode grid  -> mode_grid/
try
    g.chordFrac = linspace(0, 1.5, 20);   g.spanFrac = linspace(0, 1.5, 20);   % COARSE
    g.tspan = [0 12];   g.relTol = 1e-5;   g.absTol = 1e-7;
    g.q0 = axisAngleToQuat([0;0;1], pi/6); g.omega0 = [0;0;0];
    g.metricOpts.windowStartFrac = 0.5;    g.metricOpts.convergeTol = 0.20;
    runCoarseGrid(g, base, baseBsp, gridDir);
    fprintf('  mode grid done\n');
catch ME
    fprintf(2,'  MODE GRID FAILED: %s\n', ME.message);
end

%% 3. CoM movement suite  -> com_movement/
try
    cm = base;
    cm.tspan_pad = 0;  cm.dwellTime = 1.0;  cm.moveTime = 0.8;  cm.dt = 0.02;
    cm.odeOpts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    cm.animFps = 30;   cm.th = modeThresh; cm.animPlaybackSpeed = 0.25;
    hc = base.chordLength/2;  hs = base.spanLength/2;
    scen(1) = struct('name','S1_spanwise_sweep', ...
        'pos',{{[0;0;0],[0;0;+1.0*hs],[0;0;0],[0;0;-1.0*hs],[0;0;0]}});
    scen(2) = struct('name','S2_chordwise_sweep', ...
        'pos',{{[0;0;0],[+1.0*hc;0;0],[0;0;0],[-1.0*hc;0;0],[0;0;0]}});
    scen(3) = struct('name','S3_glide_spin_dive', ...
        'pos',{{[0.8*hc;0;0],[0;0;0],[0;0;1.2*hs],[1.5*hc;0;0]}});
    for si = 1:numel(scen)
        runComScenarioBaseline(scen(si), cm, baseBsp, comDir);
    end
    fprintf('  CoM movement suite done\n');
catch ME
    fprintf(2,'  COM SUITE FAILED: %s\n', ME.message);
end

%% 4. Test suite  -> test_suite/  (reuses seedTestCases / runOneSeedCase)
try
    ts = base;
    ts.nutPos = [0;0;0];  ts.tSamples = 0;
    ts.tspan = [0 5];  ts.odeRelTol = 1e-6;  ts.odeAbsTol = 1e-8;
    ts.enableSpanForce = true;
    ts.nIncr = 5;  ts.nutMassFrac = 0.20;  ts.nutPosMaxFrac = 1.20;
    ts.tiltMaxDeg = 45;  ts.yawSpinMin = 1;  ts.yawSpinMax = 5;
    ts.asymFactor = 0.5;  ts.stripCounts = [1 5 10 20];
    ts.groups = struct('nutMass',true,'nutChord',true,'nutSpan',true,'nutDiag',true, ...
                       'pitch',true,'roll',true,'yawSpin',true,'asymmetry',true,'stripConv',true);
    ts.savePng = true;  ts.saveOverlay = true;  ts.outputRoot = tsDir;
    runTestSuiteBaseline(ts, tsDir);
    fprintf('  test suite done\n');
catch ME
    fprintf(2,'  TEST SUITE FAILED: %s\n', ME.message);
end

fprintf('\nBASELINE SNAPSHOT COMPLETE: %s\n', snapDir);


% =========================================================================
% LOCAL: coarse mode grid (parfor) -> save modeIdx + metrics + map + ascii
% =========================================================================
function runCoarseGrid(g, base, baseBsp, outDir)
    c = base.chordLength;  S = base.spanLength;
    cfg = struct('spanLength',base.spanLength,'chordLength',base.chordLength, ...
                 'thickness',base.thickness,'bulkDensity',base.bulkDensity, ...
                 'numStrips',base.numStrips,'tSamples',0,'nutMass',base.nutMass, ...
                 'rhoFluid',base.rhoFluid,'g',base.g);
    th = defaultModeThresholds();  [modeList,modeCol] = seedModeColors();
    odeOpts = odeset('RelTol',g.relTol,'AbsTol',g.absTol);
    Nx = numel(g.chordFrac);  Nz = numel(g.spanFrac);  total = Nx*Nz;
    mIdx=zeros(total,1); dSpd=nan(total,1); vSpin=nan(total,1); tiltSd=nan(total,1);
    cone=nan(total,1); glide=nan(total,1); helixR=nan(total,1); conv=false(total,1);
    try if isempty(gcp('nocreate')); parpool('local'); end; catch; end
    parfor k=1:total
        iz=floor((k-1)/Nx)+1; ix=mod(k-1,Nx)+1;
        bsp=baseBsp; bsp.tSamples=0; bsp.nutPos_t=[g.chordFrac(ix)*c;0;g.spanFrac(iz)*S]; bsp.nutMass_t=cfg.nutMass;
        sp=buildSeedParams(bsp,cfg); x0=[zeros(3,1);g.q0;zeros(3,1);g.omega0];
        lbl='failed'; mm=[];
        try
            rhs=seedRHS(sp);
            [t,x]=ode45(@(tt,xx)rhs(tt,xx,sp),g.tspan,x0,odeOpts);
            if any(~isfinite(x(:))); lbl='chaotic'; else; mm=computeTrajectoryMetrics(t,x,g.metricOpts); lbl=classifyFlightMode(mm,th); end
        catch; lbl='failed'; end
        j=find(strcmp(lbl,modeList),1); if isempty(j); j=numel(modeList); end
        mIdx(k)=j;
        if ~isempty(mm); dSpd(k)=mm.descentSpeed; vSpin(k)=mm.verticalSpinMag; tiltSd(k)=mm.tiltStd; cone(k)=mm.coneAngleDeg; glide(k)=mm.glideRatio; helixR(k)=mm.helixRadius; conv(k)=mm.converged; end
    end
    rs=@(v) reshape(v,Nx,Nz).';
    modeIdx=rs(mIdx); metrics.descentSpeed=rs(dSpd); metrics.verticalSpinMag=rs(vSpin);
    metrics.tiltStd=rs(tiltSd); metrics.coneAngleDeg=rs(cone); metrics.glideRatio=rs(glide);
    metrics.helixRadius=rs(helixR); metrics.converged=rs(double(conv));
    chordFrac=g.chordFrac; spanFrac=g.spanFrac; q0=g.q0; omega0=g.omega0;
    save(fullfile(outDir,'mode_grid.mat'),'modeIdx','metrics','chordFrac','spanFrac','cfg','q0','omega0','g');
    % ascii map
    modeCode='gdpfstacux';
    fid=fopen(fullfile(outDir,'mode_grid_ascii.txt'),'w');
    for iz=Nz:-1:1; fprintf(fid,'span %5.2f | %s\n', spanFrac(iz), modeCode(modeIdx(iz,:))); end
    fclose(fid);
    % map png
    sc=base.spanLength/base.chordLength; xSpan=spanFrac*sc; yChord=chordFrac;
    MI=modeIdx.'; img=cat(3, reshape(modeCol(MI,1),Nx,Nz), reshape(modeCol(MI,2),Nx,Nz), reshape(modeCol(MI,3),Nx,Nz));
    f=figure('Color','w','Visible','off','Position',[80 80 900 700]);
    image(xSpan,yChord,img); set(gca,'YDir','normal'); axis image; hold on;
    present=unique(modeIdx(:)).';
    for j=present; plot(nan,nan,'s','MarkerFaceColor',modeCol(j,:),'MarkerEdgeColor','k','MarkerSize',9,'DisplayName',modeList{j}); end
    legend('Location','eastoutside'); xlabel('spanwise nut offset (x chord)'); ylabel('chordwise nut offset (x chord)');
    title(sprintf('Coarse mode grid  %dx%d', Nz, Nx));
    exportgraphics(f, fullfile(outDir,'mode_grid.png'),'Resolution',130); close(f);
end


% =========================================================================
% LOCAL: one CoM-movement scenario -> trajectory + per-dwell modes + animation
% =========================================================================
function runComScenarioBaseline(scen, cm, baseBsp, outDir)
    posList = [scen.pos{:}];
    dwell = cm.dwellTime*ones(1,size(posList,2));
    [tD,nutPos,eventT,dwellInt] = buildComPathV(posList, dwell, cm.moveTime, cm.dt);
    cfg = struct('rhoFluid',cm.rhoFluid,'g',cm.g);   % physics defaults apply
    bsp=baseBsp; bsp.tSamples=tD; bsp.nutPos_t=nutPos; bsp.nutMass_t=cm.nutMass*ones(size(tD));
    sp=buildSeedParams(bsp,cfg);
    x0=[zeros(3,1);[1;0;0;0];zeros(3,1);zeros(3,1)];
    rhs=seedRHS(sp);
    [t,x]=ode45(@(tt,xx)rhs(tt,xx,sp),[tD(1) tD(end)],x0,cm.odeOpts);
    % per-dwell settled classification
    dwellModes = strings(size(dwellInt,1),1);
    for i=1:size(dwellInt,1)
        sel = t>=dwellInt(i,1)+0.45*(dwellInt(i,2)-dwellInt(i,1)) & t<=dwellInt(i,2);
        if nnz(sel)>=8; m=computeTrajectoryMetrics(t(sel),x(sel,:),struct('windowStartFrac',0)); mf=m; mf.converged=true; dwellModes(i)=classifyFlightMode(mf,cm.th); else; dwellModes(i)="(short)"; end
    end
    save(fullfile(outDir,[scen.name '.mat']),'t','x','tD','nutPos','eventT','dwellInt','dwellModes');
    % animation + snapshot
    vfile=fullfile(outDir,[scen.name '.mp4']);
    animateModeTrajectory(t,x,struct('videoFile',vfile,'fps',cm.animFps,'title',scen.name, ...
        'eventTimes',eventT,'playbackSpeed',cm.animPlaybackSpeed,'seedParams',sp,'showSeedVels',false));
    exportgraphics(gcf, fullfile(outDir,[scen.name '_snapshot.png']),'Resolution',120);
    close(gcf);
end


% =========================================================================
% LOCAL: sweep test suite (reuses seedTestCases / runOneSeedCase)
% =========================================================================
function runTestSuiteBaseline(cfg, runDir)
    [cases, baselineBsp] = seedTestCases(cfg);
    baselineSeedParams = buildSeedParams(baselineBsp, cfg);
    results = [];
    for k=1:numel(cases)
        r = runOneSeedCase(cases(k), baselineSeedParams, baselineBsp, cfg);
        caseDir = fullfile(runDir, r.group); if ~exist(caseDir,'dir'); mkdir(caseDir); end
        result = r; save(fullfile(caseDir,[r.label '.mat']),'result');
        if ~isempty(r.x)
            fig = plotSeedCaseTrajectory(r,'off'); savefig(fig, fullfile(caseDir,[r.label '.fig']));
            exportgraphics(fig, fullfile(caseDir,[r.label '.png']),'Resolution',120); close(fig);
        end
        results=[results; r]; %#ok<AGROW>
    end
    groupNames=unique({results.group},'stable');
    for gi=1:numel(groupNames)
        gRes=results(strcmp({results.group},groupNames{gi}));
        fig=plotGroupOverlay(gRes,groupNames{gi},'off');
        savefig(fig, fullfile(runDir,groupNames{gi},'_overlay.fig'));
        exportgraphics(fig, fullfile(runDir,groupNames{gi},'_overlay.png'),'Resolution',120); close(fig);
    end
    caseSummary=struct('group',{results.group},'label',{results.label},'status',{results.status}, ...
                       'sweepVal',{results.sweepVal},'descended',{results.descended});
    save(fullfile(runDir,'manifest.mat'),'cfg','caseSummary');
    fid=fopen(fullfile(runDir,'summary.txt'),'w');
    fprintf(fid,'%-11s %-24s %-10s %-9s\n','group','label','status','descended');
    for k=1:numel(results); fprintf(fid,'%-11s %-24s %-10s %-9d\n',results(k).group,results(k).label,results(k).status,results(k).descended); end
    fclose(fid);
end


% =========================================================================
% LOCAL: waypoint list (+ dwells) -> dense pchip nut path + dwell intervals
% =========================================================================
function [tDense,nutPos,eventT,dwellInt] = buildComPathV(posList,dwellList,moveTime,dt)
    K=size(posList,2); wpT=0; wpP=posList(:,1); tcur=0; eventT=zeros(1,K); dwellInt=zeros(K,2);
    for i=1:K
        eventT(i)=tcur; dwellInt(i,1)=tcur;
        if dwellList(i)>0; tcur=tcur+dwellList(i); wpT(end+1)=tcur; wpP(:,end+1)=posList(:,i); end %#ok<AGROW>
        dwellInt(i,2)=tcur;
        if i<K; tcur=tcur+moveTime; wpT(end+1)=tcur; wpP(:,end+1)=posList(:,i+1); end %#ok<AGROW>
    end
    tDense=0:dt:wpT(end); if tDense(end)<wpT(end); tDense(end+1)=wpT(end); end
    nutPos=zeros(3,numel(tDense));
    for a=1:3; nutPos(a,:)=interp1(wpT,wpP(a,:),tDense,'pchip'); end
end


% =========================================================================
% LOCAL: small utilities
% =========================================================================
function s = gitField(repoRoot, args)
    [st,out] = system(sprintf('git -C "%s" %s', repoRoot, args));
    if st==0; s=strtrim(out); else; s='nogit'; end
end
function writeManifestTxt(path, m)
    fid=fopen(path,'w');
    fprintf(fid,'Model baseline snapshot\n=======================\n');
    fprintf(fid,'timestamp : %s\n', m.stamp);
    fprintf(fid,'git commit: %s%s\n', m.gitHash, ternary(m.gitDirty,'  (WORKING TREE DIRTY)',''));
    fprintf(fid,'matlab    : %s\n\n', m.matlab);
    fprintf(fid,'Base seed:\n'); disp_struct(fid, m.baseSeed);
    fprintf(fid,'\nPhysics switches:\n'); disp_struct(fid, m.switches);
    fprintf(fid,'\nAero coeffs (at alpha=0; law constants pinned by git hash):\n'); disp_struct(fid, m.aeroCoeffs);
    fprintf(fid,'\nMode-classifier thresholds:\n'); disp_struct(fid, m.modeThresholds);
    fclose(fid);
end
function disp_struct(fid, s)
    f=fieldnames(s);
    for i=1:numel(f)
        v=s.(f{i});
        if isnumeric(v)&&isscalar(v); fprintf(fid,'  %-28s = %g\n', f{i}, v);
        elseif islogical(v)&&isscalar(v); fprintf(fid,'  %-28s = %d\n', f{i}, v);
        elseif isnumeric(v); fprintf(fid,'  %-28s = [%s]\n', f{i}, num2str(v(:).'));
        else; fprintf(fid,'  %-28s = <%s>\n', f{i}, class(v)); end
    end
end
function out = ternary(cond, a, b); if cond; out=a; else; out=b; end; end
