function res = testRefactorInvariance(mode, opts)
% TESTREFACTORINVARIANCE  Guard for code cleanups: the model must not change.
%
% A refactor (moving the builder, archiving the planar model, retiring the old
% classifier, ...) should leave the PHYSICS and the LABELS exactly as they were.
% This captures a "golden" reference from the code as it stands, and later
% compares the current code against it, bit for bit.
%
%   testRefactorInvariance('capture')   write model_test_results/fixtures/golden_<git>.mat
%   testRefactorInvariance('compare')   compare against the newest golden_*.mat
%   testRefactorInvariance('compare', struct('fixture', file))   ... against a given one
%
% WHAT IS PINNED
%   1. The built seedParams struct of a fixed set of seeds (flat; twist 5/20 deg;
%      bowl and one-sided curvature; a moving nut; the Hou et al. plate at its
%      five anchors; LEV colocated/forward; the 3D added-mass tensor). Compared
%      field by field: numbers with isequaln, function handles by their text,
%      polyshapes by their vertices.
%   2. seed6DOFODE3D at 50 fixed random states per seed: xdot, the aero force and
%      the aero torque. Must be BIT-identical (max abs difference exactly 0).
%   3. Twelve 4 s trajectories: the final state (bit-identical, so the whole
%      integration is pinned), the six-part label and parts, the paper-vocabulary
%      label, and descentSpeed / absSpanTiltDeg / revolutionRate (1e-12 rel).
%
% A difference is reported with the path of the field or the seed/run it is in.
% Fields a cleanup removes on purpose are listed in opts.allowRemoved (cellstr of
% dotted paths below seedParams, e.g. {'enableSpanForce'}), and are then only
% reported, not failed.
%
% OUTPUT
%   res : struct array (.name, .pass, .detail), one entry per check.

    if nargin < 1 || isempty(mode); mode = 'compare'; end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'allowRemoved'); opts.allowRemoved = {}; end

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  P = seedPaths(here);
    fixDir = fullfile(P.results, 'fixtures');

    cur = computeReference();

    switch lower(mode)
        case 'capture'
            if ~exist(fixDir, 'dir'); mkdir(fixDir); end
            [~, gh] = system(sprintf('git -C "%s" rev-parse --short HEAD', P.root));
            [~, gd] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', P.root));
            tag = strtrim(gh);  if ~isempty(strtrim(gd)); tag = [tag '-dirty']; end
            file = fullfile(fixDir, sprintf('golden_%s_%s.mat', ...
                            datestr(now, 'yyyy-mm-dd_HHMMSS'), tag)); %#ok<TNOW1,DATST>
            golden = cur; %#ok<NASGU>
            save(file, 'golden');
            fprintf('captured %d seeds, %d trajectories -> %s\n', ...
                    numel(cur.seeds), numel(cur.runs), file);
            res = struct('name', 'capture', 'pass', true, 'detail', file);
            return
        case 'compare'
            if isfield(opts, 'fixture')
                file = opts.fixture;
            else
                d = dir(fullfile(fixDir, 'golden_*.mat'));
                if isempty(d)
                    error('testRefactorInvariance:noFixture', ...
                          'No golden_*.mat in %s: run testRefactorInvariance(''capture'') first.', fixDir);
                end
                [~, o] = sort({d.name});  d = d(o);
                file = fullfile(d(end).folder, d(end).name);
            end
            S = load(file);  gold = S.golden;
            fprintf('comparing against %s\n', file);
        otherwise
            error('testRefactorInvariance:mode', 'mode must be ''capture'' or ''compare''.');
    end

    res = struct('name', {}, 'pass', {}, 'detail', {});

    % --- 1. built seedParams ------------------------------------------------
    nameG = {gold.seeds.name};
    for k = 1:numel(cur.seeds)
        g = find(strcmp(nameG, cur.seeds(k).name), 1);
        if isempty(g)
            res(end+1) = item(['seed ' cur.seeds(k).name], false, 'not in the fixture'); %#ok<AGROW>
            continue
        end
        [diffs, removed] = structDiff(gold.seeds(g).sp, cur.seeds(k).sp, 'sp', opts.allowRemoved);
        det = '';
        if ~isempty(removed); det = ['removed as allowed: ' strjoin(removed, ', ')]; end
        if ~isempty(diffs); det = strjoin([diffs(1:min(end,6)), {det}], ' | '); end
        res(end+1) = item(['seedParams ' cur.seeds(k).name], isempty(diffs), det); %#ok<AGROW>

        dx = max(abs(gold.seeds(g).xdot(:) - cur.seeds(k).xdot(:)));
        dF = max(abs(gold.seeds(g).F(:)    - cur.seeds(k).F(:)));
        dT = max(abs(gold.seeds(g).tau(:)  - cur.seeds(k).tau(:)));
        res(end+1) = item(['RHS ' cur.seeds(k).name], dx == 0 && dF == 0 && dT == 0, ...
                          sprintf('max|dxdot| %.1e  |dF| %.1e  |dtau| %.1e', dx, dF, dT)); %#ok<AGROW>
    end

    % --- 2. trajectories, labels, metrics ------------------------------------
    nameR = {gold.runs.name};
    for k = 1:numel(cur.runs)
        c = cur.runs(k);   g = find(strcmp(nameR, c.name), 1);
        if isempty(g)
            res(end+1) = item(['run ' c.name], false, 'not in the fixture'); %#ok<AGROW>
            continue
        end
        o = gold.runs(g);
        dEnd = max(abs(o.xEnd(:) - c.xEnd(:)));
        relM = max(abs([o.descentSpeed o.absSpanTiltDeg o.revolutionRate] - ...
                       [c.descentSpeed c.absSpanTiltDeg c.revolutionRate]) ./ ...
                   max(abs([o.descentSpeed o.absSpanTiltDeg o.revolutionRate]), eps));
        ok = dEnd == 0 && strcmp(o.label, c.label) && strcmp(o.parts, c.parts) && ...
             strcmp(o.paperLabel, c.paperLabel) && ~(relM > 1e-12);
        res(end+1) = item(['run ' c.name], ok, sprintf( ...
            '%s [%s] / %s  vs gold %s [%s] / %s;  |dx_end| %.1e, metrics rel %.1e', ...
            c.label, c.parts, c.paperLabel, o.label, o.parts, o.paperLabel, dEnd, relM)); %#ok<AGROW>
    end

    % --- report ---------------------------------------------------------------
    nPass = nnz([res.pass]);
    for k = 1:numel(res)
        if res(k).pass; tag = 'PASS'; else; tag = 'FAIL'; end
        fprintf('  %-34s %s  %s\n', res(k).name, tag, res(k).detail);
    end
    if nPass == numel(res)
        fprintf('PASS: %d/%d -- the model is unchanged by the refactor.\n', nPass, numel(res));
    else
        fprintf(2, 'FAIL: %d/%d checks differ from the golden reference.\n', ...
                numel(res) - nPass, numel(res));
    end
end


% =========================================================================
% The reference itself: a fixed seed set, RHS samples, trajectories
% =========================================================================
function ref = computeReference()
    seeds = seedSet();
    ref.seeds = struct('name', {}, 'sp', {}, 'xdot', {}, 'F', {}, 'tau', {});
    for k = 1:numel(seeds)
        sp = buildSeedParams(seeds(k).bsp, seeds(k).cfg);
        rng(7);                                    % same 50 states for every seed
        N = 50;   tMax = max(sp.massParams.tSamples(end), 0);
        xd = zeros(13, N);  F = zeros(3, N);  T = zeros(3, N);
        for j = 1:N
            q = randn(4,1);  q = q / norm(q);
            x = [randn(3,1); q; 2*randn(3,1); 10*randn(3,1)];
            t = tMax * rand();
            [xd(:,j), im] = seed6DOFODE3D(t, x, sp);
            F(:,j) = im.F_aero_body;   T(:,j) = im.tau_body;
        end
        ref.seeds(k) = struct('name', seeds(k).name, 'sp', sp, 'xdot', xd, 'F', F, 'tau', T);
    end

    th = defaultSeedModeThresholds();
    ref.runs = struct('name', {}, 'xEnd', {}, 'label', {}, 'parts', {}, 'paperLabel', {}, ...
                      'descentSpeed', {}, 'absSpanTiltDeg', {}, 'revolutionRate', {});
    for k = find([seeds.run])
        s  = seeds(k);
        sp = buildSeedParams(s.bsp, s.cfg);
        [t, x] = ode45(@(tt,xx) seed6DOFODE3D(tt, xx, sp), [0 4], ...
                       [zeros(3,1); s.q0; zeros(3,1); s.w0], odeset('RelTol',1e-6,'AbsTol',1e-8));
        mo = s.metricOpts;   mo.flipHysteresisTurns = th.flipHysteresisTurns;
        sm = computeSeedModeMetrics(t, x, mo);
        [lbl, p] = classifySeedMode(sm, th);
        ref.runs(end+1) = struct('name', s.name, 'xEnd', x(end,:), 'label', lbl, ...
            'parts', sprintf('%s/%s/%s/%s/%s/%d', p.flip, p.revolution, p.attitude, ...
                             p.path, p.regularity, p.physical), ...
            'paperLabel', classifyPaperMode(sm), 'descentSpeed', sm.descentSpeed, ...
            'absSpanTiltDeg', sm.absSpanTiltDeg, 'revolutionRate', sm.revolutionRate);
    end
end

function seeds = seedSet()
% The fixed seed set. Our working seed (the suite's), plus the Hou et al. plate.
    S = 0.050;  c = 0.015;  th = 0.002;  hs = S/2;  hc = c/2;
    base.seedShape     = polyshape([-hs, hs, hs, -hs], [-hc, -hc, hc, hc]);
    base.seedDensity   = 65 * th;
    base.seedThickness = th;
    base.numStrips     = 10;
    base.tSamples      = 0;
    base.nutMass_t     = 75e-6;
    base.nutPos_t      = [0;0;0];
    cfg  = struct('rhoFluid', 1.225, 'g', 9.81, 'shapeModel', 'shape3d');
    q0   = axisAngleToQuat([0;0;1], pi/6);
    mo   = struct('windowStartFrac', 0.5, 'convergeTol', 0.20, 'spanLength', S);

    seeds = struct('name', {}, 'bsp', {}, 'cfg', {}, 'run', {}, 'q0', {}, 'w0', {}, 'metricOpts', {});
    add = @(seeds, nm, b, cf, run, q) [seeds, struct('name', nm, 'bsp', b, 'cfg', cf, ...
                                       'run', run, 'q0', q, 'w0', [0;0;0], 'metricOpts', mo)];

    seeds = add(seeds, 'flat_centre', base, cfg, true, q0);
    b = base;  b.nutPos_t = [0.5*hc; 0; 0.6*hs];
    seeds = add(seeds, 'flat_offset', b, cfg, true, q0);
    for tip = [5 20]
        b = base;  kt = deg2rad(tip)/hs;  b.twist = @(z) kt .* z;
        seeds = add(seeds, sprintf('twist_%02d', tip), b, cfg, true, [1;0;0;0]);
    end
    kc = deg2rad(20)/hs;
    b = base;  b.curvature = @(s) kc .* s;
    seeds = add(seeds, 'bowl_20', b, cfg, true, q0);
    b = base;  b.curvature = @(s) (s > 0) .* kc .* s;
    seeds = add(seeds, 'oneSide_20', b, cfg, true, q0);
    cf = cfg;  cf.enableAddedMass3D = true;
    b = base;  b.curvature = @(s) kc .* s;
    seeds = add(seeds, 'bowl_20_addedMass3D', b, cf, false, q0);

    % moving nut (the suite's stage-D scenario 3)
    [tD, nutPos] = buildComPath([[0.8*hc;0;0], [0;0;0], [0;0;1.2*hs], [1.5*hc;0;0]], ...
                                [1 1 1 1], 0.3, 0.01);
    b = base;  b.tSamples = tD;  b.nutPos_t = nutPos;  b.nutMass_t = 75e-6 * ones(size(tD));
    seeds = add(seeds, 'comMove', b, cfg, true, q0);

    % LEV, both application points
    b = base;  b.nutPos_t = [0.5*hc; 0; 0.6*hs];
    cf = cfg;  cf.enableLEV = true;  cf.levApplicationPoint = 'colocated';
    seeds = add(seeds, 'lev_colocated', b, cf, true, q0);
    cf.levApplicationPoint = 'forward';
    seeds = add(seeds, 'lev_forward', b, cf, true, q0);
    b = base;  kt = deg2rad(20)/hs;  b.twist = @(z) kt .* z;
    seeds = add(seeds, 'twist_20_lev_forward', b, cf, false, [1;0;0;0]);

    % the Hou et al. plate at its five published anchors
    [pcfg, pbsp, toNutPos, pinfo] = paperSeedConfig();
    pmo = struct('windowStartFrac', 0.5, 'convergeTol', 0.20, 'refLength', pinfo.L, ...
                 'sigma', pinfo.sigma, 'g', pcfg.g, 'rhoFluid', pcfg.rhoFluid, ...
                 'spanLength', pinfo.L);
    anch = {'AR', 0.39, 0.17; 'CST', 0.06, 0.04; 'SST', 0.10, 0.10; 'CH', 0.18, 0.08; 'FA', 0.35, 0.08};
    for i = 1:size(anch, 1)
        b = pbsp;  b.tSamples = 0;
        b.nutPos_t = toNutPos(anch{i,2}, anch{i,3});  b.nutMass_t = pcfg.nutMass;
        s = struct('name', ['paper_' anch{i,1}], 'bsp', b, 'cfg', pcfg, 'run', i <= 3, ...
                   'q0', pcfg.releaseQuat, 'w0', pcfg.releaseOmega, 'metricOpts', pmo);
        seeds = [seeds, s]; %#ok<AGROW>
    end
end


% =========================================================================
% Field-by-field comparison of two values
% =========================================================================
function [diffs, removed] = structDiff(a, b, path, allowRemoved)
    diffs = {};  removed = {};
    if isstruct(a) && isstruct(b)
        if ~isequal(size(a), size(b))
            diffs{end+1} = sprintf('%s: struct size %s vs %s', path, mat2str(size(a)), mat2str(size(b)));
            return
        end
        fa = fieldnames(a);  fb = fieldnames(b);
        for f = setdiff(fa, fb).'
            p = [path '.' f{1}];
            if any(strcmp(strrep(p, 'sp.', ''), allowRemoved)); removed{end+1} = p; %#ok<AGROW>
            else; diffs{end+1} = [p ': missing now']; end %#ok<AGROW>
        end
        for f = setdiff(fb, fa).'
            diffs{end+1} = [path '.' f{1} ': new field']; %#ok<AGROW>
        end
        for f = intersect(fa, fb).'
            for e = 1:numel(a)
                [d, r] = structDiff(a(e).(f{1}), b(e).(f{1}), [path '.' f{1}], allowRemoved);
                diffs = [diffs, d];  removed = [removed, r]; %#ok<AGROW>
            end
        end
    elseif isa(a, 'function_handle') && isa(b, 'function_handle')
        if ~strcmp(func2str(a), func2str(b)); diffs{end+1} = [path ': handle text differs']; end
    elseif isa(a, 'polyshape') && isa(b, 'polyshape')
        if ~isequaln(a.Vertices, b.Vertices); diffs{end+1} = [path ': polyshape differs']; end
    elseif iscell(a) && iscell(b) && isequal(size(a), size(b))
        for e = 1:numel(a)
            [d, r] = structDiff(a{e}, b{e}, sprintf('%s{%d}', path, e), allowRemoved);
            diffs = [diffs, d];  removed = [removed, r]; %#ok<AGROW>
        end
    elseif ~isequaln(a, b)
        if isnumeric(a) && isnumeric(b) && isequal(size(a), size(b))
            diffs{end+1} = sprintf('%s: max|diff| %.2e', path, max(abs(a(:) - b(:))));
        else
            diffs{end+1} = [path ': differs (' class(a) ' vs ' class(b) ')'];
        end
    end
end

function r = item(name, pass, detail)
    r = struct('name', name, 'pass', logical(pass), 'detail', detail);
end
