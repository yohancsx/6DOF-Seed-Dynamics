function P = seedPaths(root)
% SEEDPATHS  Put the model code on the MATLAB path and return the project's folders.
%
% Every script calls this instead of hardcoding a machine-specific path. The repo
% root is found from THIS file's location, so the code runs from any clone, on any
% machine, whatever MATLAB's current folder is.
%
% A script finds this file by walking up from its own folder (it cannot call
% seedPaths before seedPaths is on the path). The three lines every script uses:
%
%   here = fileparts(mfilename('fullpath'));  if isempty(here); here = pwd; end
%   while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
%   addpath(here);  P = seedPaths(here);
%
% (mfilename is empty when code is run with Evaluate Selection, hence the pwd
% fallback: run from inside the repo in that case.)
%
% FOLDERS ADDED TO THE PATH
%   physics/ (+ helpers/, aero/, mass/), physics3d/, visualization/, testing/helpers/
%   Scripts that need a suite-specific folder (e.g. testing/shape3d/ for
%   reportPaperModeGrid) add it themselves from P.root.
%
% WHERE OUTPUTS GO (the project rule)
%   P.results : model_test_results/ in the repo (git-ignored). Everything a test
%               suite or regression produces: runs that are regenerated whenever
%               the model changes.
%   P.misc    : "Misc Outputs" beside the repo. One-off studies and tooling
%               (classifier review sets, the mode gallery, LEV exploration, the
%               single-drop harnesses).
%   The separate "Outputs" folder beside the repo is reserved for material made for
%   the paper or a presentation. No test or study writes there.
%
% INPUT
%   root : (optional) the repo root. Default: the folder holding this file.
% OUTPUT
%   P    : struct with .root, .results, .misc, .seedData (the measured-seed shape
%          file, outside the repo; may not exist on every machine).

    if nargin < 1 || isempty(root)
        root = fileparts(mfilename('fullpath'));
    end
    root = char(root);
    if ~isfile(fullfile(root, 'seedPaths.m'))
        error('seedPaths:badRoot', ...
              'No seedPaths.m in "%s": run the script from inside the 6DOF Seed Dynamics repo.', root);
    end

    addpath(fullfile(root, 'physics'), fullfile(root, 'physics', 'helpers'), ...
            fullfile(root, 'physics', 'aero'), fullfile(root, 'physics', 'mass'), ...
            fullfile(root, 'physics3d'), fullfile(root, 'visualization'), ...
            fullfile(root, 'testing', 'helpers'));

    parent     = fileparts(root);
    P.root     = root;
    P.results  = fullfile(root, 'model_test_results');
    P.misc     = fullfile(parent, 'Misc Outputs');
    P.seedData = fullfile(parent, 'Real Seed Validation Data', 'seedResults.mat');
end
