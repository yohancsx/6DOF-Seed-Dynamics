function seeds = loadRealSeeds(file, opts)
% LOADREALSEEDS  The measured seeds of seedMassModel.mat, ready for buildSeedParams.
%
%   seeds = loadRealSeeds()                       all seeds, P.seedMassModel, defaults
%   seeds = loadRealSeeds(file, struct('seeds', [5 24], 'calibration', 'nutOnly'))
%
% Each seed becomes base seed params carrying an EXPLICIT strip geometry and mass
% model (bsp.stripGeom; see realSeedToBsp), which setupSeedShape3D takes as given:
%   sp = buildSeedParams(seeds(k).bsp, cfg);
% The strips are the CT-fitted plates; the mass, CoM and inertia are the fitted
% whole-seed values (nut + wing). Nothing is tuned to flight data.
%
% Recommended cfg for these seeds: cfg.enableAddedMass3D = true. Their plates
% carry real dihedral (tips to ~26 deg), and the default flat-plate added-mass form
% assumes every strip normal is body y.
%
% INPUTS
%   file : (optional) the .mat; default P.seedMassModel from seedPaths.
%   opts : (optional) struct:
%            .seeds       seed numbers to load (default: all, in file order)
%            .calibration 'both' (default) | 'nutOnly'  -- see realSeedToBsp
%            .subStrips   strips per plate (default 1); cuts along the plate span
% OUTPUT
%   seeds : struct array, one per seed: .name, .seedNum, .bsp, .meta.

    if nargin < 1 || isempty(file)
        here = fileparts(mfilename('fullpath'));
        while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
        P = seedPaths(here);  file = P.seedMassModel;
    end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'calibration'); opts.calibration = 'both'; end
    if ~isfield(opts, 'subStrips');   opts.subStrips   = 1;      end
    opts.calibration = char(opts.calibration);
    if ~any(strcmp(opts.calibration, {'both', 'nutOnly'}))
        error('loadRealSeeds:badCalibration', 'opts.calibration must be ''both'' or ''nutOnly''.');
    end
    if ~(isscalar(opts.subStrips) && opts.subStrips >= 1 && opts.subStrips == round(opts.subStrips))
        error('loadRealSeeds:badSubStrips', 'opts.subStrips must be a positive integer.');
    end
    if ~isfile(file)
        error('loadRealSeeds:noFile', 'Seed mass model not found: %s', file);
    end

    S = load(file, 'out', 'TP');
    nums = [S.out.seedNum];
    if isfield(opts, 'seeds') && ~isempty(opts.seeds)
        missing = setdiff(opts.seeds, nums);
        if ~isempty(missing)
            error('loadRealSeeds:unknownSeed', 'Seed(s) %s not in %s.', mat2str(missing), file);
        end
        idx = arrayfun(@(s) find(nums == s, 1), opts.seeds);
    else
        idx = 1:numel(S.out);
    end

    seeds = struct('name', {}, 'seedNum', {}, 'bsp', {}, 'meta', {});
    for k = idx
        o = S.out(k);
        parts = S.TP(S.TP.seed == string(o.name), :);
        [bsp, meta] = realSeedToBsp(o, parts, opts);
        seeds(end+1) = struct('name', o.name, 'seedNum', o.seedNum, 'bsp', bsp, 'meta', meta); %#ok<AGROW>
    end
end
