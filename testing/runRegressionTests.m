function [T, summary] = runRegressionTests(opts)
% RUNREGRESSIONTESTS  Run every pass/fail regression test, in one go.
%
% The fast checks that should hold after ANY change (about a minute in all):
%   testMinimalImpEquivalence   the 3D model in 2D == Olivia's minimal_imp (primary
%                               flat test; its CoM-offset cases WARN until the
%                               added-mass coupling is in, see the README roadmap)
%   testShape3DFlatEquivalence  the 3D model == the archived planar model's core
%   testCurvedMassProperties    curved-seed mass and inertia
%   testNewPhysicsTerms         added-mass rate + Munk moment, edge drag, LEV
%   testPaperMetrics            paper metrics and both classifiers on synthetic motion
%   testActivePathClean         no active code calls archived code
%
% NOT included: testRefactorInvariance, which pins the model bit for bit against a
% captured reference. It is for refactors only -- an intended physics change
% fails it by design. Run it by hand around a cleanup.
%
% A test that throws is reported as one failed check, so a broken test never stops
% the rest. Called by stage 0 of runSeedTestSuite3D, which only WARNS on failure.
%
% INPUT   opts.tests : (optional) cellstr subset of the test names above.
% OUTPUT  T       : table, one row per check: test, check, pass, warn, detail.
%         summary : struct with .nPass, .nTotal, .nWarn, .failed (string list),
%                   .lines (string column for a summary file).

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  P = seedPaths(here);
    addpath(fullfile(P.root, 'testing', 'shape3d'));

    tests = {'testMinimalImpEquivalence', 'testShape3DFlatEquivalence', ...
             'testCurvedMassProperties', 'testNewPhysicsTerms', 'testPaperMetrics', ...
             'testActivePathClean'};
    if nargin >= 1 && isfield(opts, 'tests'); tests = cellstr(opts.tests); end

    rows = cell(0, 5);
    for k = 1:numel(tests)
        fprintf('\n######## %s\n', tests{k});
        t0 = tic;
        try
            r = feval(tests{k});
            for j = 1:numel(r)
                rows(end+1, :) = {string(tests{k}), string(r(j).name), logical(r(j).pass), ...
                                  startsWith(string(r(j).detail), "WARN"), string(r(j).detail)}; %#ok<AGROW>
            end
        catch ME
            rows(end+1, :) = {string(tests{k}), "(test threw an error)", false, false, ...
                              string(ME.message)}; %#ok<AGROW>
            fprintf(2, 'ERROR in %s: %s\n', tests{k}, ME.message);
        end
        fprintf('(%s: %.0f s)\n', tests{k}, toc(t0));
    end
    T = cell2table(rows, 'VariableNames', {'test', 'check', 'pass', 'warn', 'detail'});

    summary.nTotal = height(T);
    summary.nPass  = nnz(T.pass);
    summary.nWarn  = nnz(T.warn & T.pass);
    bad = T(~T.pass, :);
    summary.failed = bad.test + ": " + bad.check;
    L = strings(0, 1);
    L(end+1) = sprintf('REGRESSION TESTS: %d/%d pass, %d warn', summary.nPass, summary.nTotal, summary.nWarn);
    for name = unique(T.test, 'stable').'
        R = T(T.test == name, :);
        L(end+1) = sprintf('  %-28s %d/%d pass', name, nnz(R.pass), height(R)); %#ok<AGROW>
        for k = find(~R.pass | R.warn).'
            if ~R.pass(k); tag = "FAIL"; else; tag = "WARN"; end
            L(end+1) = sprintf('      %s  %s -- %s', tag, R.check(k), R.detail(k)); %#ok<AGROW>
        end
    end
    summary.lines = L(:);
    fprintf('\n%s\n', strjoin(summary.lines, newline));
end
