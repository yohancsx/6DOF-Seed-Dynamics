function compareBaselineModels(snapDir, refSnapDir)
%COMPAREBASELINEMODELS  Diff the planar and shape3d stages of a baseline snapshot.
%
% generateModelBaseline runs the three shared stages (mode grid, CoM movement,
% test suite) under BOTH models on the SAME flat seed, writing them to
% mode_grid/ + mode_grid_3d/ (and likewise for the others). Because the seed is
% identical, every difference between a pair of folders is the MODEL: 'planar'
% still carries Tx, Ty and the span torque; 'shape3d' is the cleaned model.
%
% This is the comparison to iterate against. Re-run it after each phase-4 change
% and watch the mode-migration table: the goal is for shape3d to recover (and
% then improve on) the planar mode structure through real physics rather than
% through the terms the Sep-2026 audit removed.
%
% USAGE
%   compareBaselineModels                      % newest snapshot, planar vs shape3d
%   compareBaselineModels(snapDir)             % that snapshot, planar vs shape3d
%   compareBaselineModels(snapDir, refSnapDir) % also regression-check the planar
%                                              % stages against an older snapshot
%
% INPUTS
%   snapDir    : snapshot folder (default: newest under model_test_results/).
%   refSnapDir : (optional) older snapshot whose PLANAR stages should be
%                byte-identical -- the frozen-model regression check.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    resRoot = fullfile(root, 'model_test_results');
    addpath(fullfile(root,'testing','helpers'));

    if nargin < 1 || isempty(snapDir)
        d = dir(fullfile(resRoot, '20*'));  d = d([d.isdir]);
        assert(~isempty(d), 'No snapshots under %s', resRoot);
        [~, newest] = max([d.datenum]);
        snapDir = fullfile(resRoot, d(newest).name);
    end
    fprintf('Snapshot: %s\n', snapDir);

    %% 1. Mode grid: planar vs shape3d
    fA = fullfile(snapDir,'mode_grid','mode_grid.mat');
    fB = fullfile(snapDir,'mode_grid_3d','mode_grid.mat');
    if ~isfile(fA) || ~isfile(fB)
        fprintf(2,'  mode grid pair missing -- was this snapshot built with dual models?\n');
        return
    end
    A = load(fA);  B = load(fB);
    [modeList, modeCol] = seedModeColors();

    changed = A.modeIdx ~= B.modeIdx;
    nTot = numel(changed);
    fprintf('\n=== MODE GRID: planar -> shape3d ===\n');
    fprintf('  cells: %d   changed: %d (%.1f%%)   unchanged: %.1f%%\n', ...
            nTot, nnz(changed), 100*nnz(changed)/nTot, 100*(1-nnz(changed)/nTot));

    % Mode census, both models
    fprintf('\n  %-14s %8s %8s %9s\n','mode','planar','shape3d','delta');
    for j = 1:numel(modeList)
        na = nnz(A.modeIdx==j);  nb = nnz(B.modeIdx==j);
        if na==0 && nb==0; continue; end
        fprintf('  %-14s %8d %8d %+9d\n', modeList{j}, na, nb, nb-na);
    end

    % Migration table: which planar mode became which
    fprintf('\n  Migration (planar row -> shape3d col), cells that CHANGED:\n');
    src = unique(A.modeIdx(changed)).';
    for j = src
        dst = B.modeIdx(changed & A.modeIdx==j);
        u = unique(dst).';
        parts = arrayfun(@(k) sprintf('%s x%d', modeList{k}, nnz(dst==k)), u, 'uni', 0);
        fprintf('    %-14s -> %s\n', modeList{j}, strjoin(parts, ',  '));
    end

    % Metric deltas over the cells that stayed the same mode (like-for-like)
    same = ~changed;
    fprintf('\n  Metric shift on cells that kept their mode (n=%d):\n', nnz(same));
    mf = {'descentSpeed','coneAngleDeg','verticalSpinMag','tiltStd'};
    fprintf('  %-18s %10s %10s %10s\n','metric','planar','shape3d','delta');
    for i = 1:numel(mf)
        if ~isfield(A.metrics,mf{i}); continue; end
        a = A.metrics.(mf{i})(same);  b = B.metrics.(mf{i})(same);
        ok = isfinite(a) & isfinite(b);
        fprintf('  %-18s %10.3f %10.3f %+10.3f\n', mf{i}, ...
                mean(a(ok)), mean(b(ok)), mean(b(ok)-a(ok)));
    end

    %% 2. Side-by-side map + difference
    Nz = size(A.modeIdx,1);  Nx = size(A.modeIdx,2);
    sc = A.cfg.spanLength / A.cfg.chordLength;
    xSpan = A.spanFrac*sc;  yChord = A.chordFrac;
    f = figure('Color','w','Visible','off','Position',[60 60 1500 520]);
    ttl = {'planar (frozen reference)','shape3d (cleaned model)'};
    for p = 1:2
        subplot(1,3,p);
        MI = ternaryIdx(p, A.modeIdx, B.modeIdx).';
        img = cat(3, reshape(modeCol(MI,1),Nx,Nz), reshape(modeCol(MI,2),Nx,Nz), ...
                     reshape(modeCol(MI,3),Nx,Nz));
        image(xSpan, yChord, img); set(gca,'YDir','normal'); axis image;
        xlabel('spanwise nut offset (\times chord)'); ylabel('chordwise nut offset (\times chord)');
        title(ttl{p});
    end
    subplot(1,3,3);
    imagesc(xSpan, yChord, double(changed)); set(gca,'YDir','normal'); axis image;
    colormap(gca, [0.93 0.93 0.93; 0.80 0.16 0.16]);  caxis([0 1]);
    xlabel('spanwise nut offset (\times chord)'); ylabel('chordwise nut offset (\times chord)');
    title(sprintf('changed mode (%.0f%% of cells)', 100*nnz(changed)/nTot));
    sgtitle('Mode grid: frozen planar model vs cleaned 3D model, same flat seed');
    outPng = fullfile(snapDir,'mode_grid_comparison.png');
    exportgraphics(f, outPng, 'Resolution',130);  close(f);
    fprintf('\n  wrote %s\n', outPng);

    %% 3. Test suite + CoM movement: do the two models agree per case?
    compareTextFile(snapDir, 'test_suite', 'summary.txt', 'TEST SUITE');
    compareComMovement(snapDir);

    %% 4. Optional planar regression against an older snapshot
    if nargin >= 2 && ~isempty(refSnapDir)
        fprintf('\n=== PLANAR REGRESSION vs %s ===\n', refSnapDir);
        for pair = {{'mode_grid','mode_grid_ascii.txt'}, {'test_suite','summary.txt'}}
            p = pair{1};
            fa = fullfile(refSnapDir, p{1}, p{2});   fb = fullfile(snapDir, p{1}, p{2});
            if isfile(fa) && isfile(fb)
                same = isequal(fileread(fa), fileread(fb));
                fprintf('  %-28s %s\n', fullfile(p{1},p{2}), ...
                        string(ternaryStr(same,'IDENTICAL','DIFFERS')));
            end
        end
    end
end


function out = ternaryIdx(p, a, b);  if p==1; out=a; else; out=b; end;  end
function out = ternaryStr(c, a, b);  if c;    out=a; else; out=b; end;  end


function compareTextFile(snapDir, stage, fname, label)
    fa = fullfile(snapDir, stage, fname);
    fb = fullfile(snapDir, [stage '_3d'], fname);
    if ~isfile(fa) || ~isfile(fb); return; end
    fprintf('\n=== %s: planar vs shape3d ===\n', label);
    if isequal(fileread(fa), fileread(fb))
        fprintf('  summaries IDENTICAL (the cut terms changed no case outcome here)\n');
    else
        A = splitlines(string(fileread(fa)));  B = splitlines(string(fileread(fb)));
        n = min(numel(A), numel(B));  d = find(A(1:n) ~= B(1:n));
        fprintf('  %d of %d summary lines differ; first few:\n', numel(d), n);
        for k = d(1:min(8,numel(d))).'
            fprintf('    planar : %s\n    shape3d: %s\n', strtrim(A(k)), strtrim(B(k)));
        end
    end
end


function compareComMovement(snapDir)
    scen = {'S1_spanwise_sweep','S2_chordwise_sweep','S3_glide_spin_dive'};
    fprintf('\n=== COM MOVEMENT: dwell-mode sequences ===\n');
    for i = 1:numel(scen)
        fa = fullfile(snapDir,'com_movement',[scen{i} '.mat']);
        fb = fullfile(snapDir,'com_movement_3d',[scen{i} '.mat']);
        if ~isfile(fa) || ~isfile(fb); continue; end
        a = load(fa,'dwellModes');  b = load(fb,'dwellModes');
        tag = ternaryStr(isequal(a.dwellModes,b.dwellModes), 'same', 'CHANGED');
        fprintf('  %-20s %-8s\n     planar : %s\n     shape3d: %s\n', scen{i}, tag, ...
                strjoin(cellstr(a.dwellModes(:).'),' -> '), ...
                strjoin(cellstr(b.dwellModes(:).'),' -> '));
    end
end
