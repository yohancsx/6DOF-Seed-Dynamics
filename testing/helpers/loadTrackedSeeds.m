function tracks = loadTrackedSeeds(file, opts)
% LOADTRACKEDSEEDS  The measured drop trajectories (seedTrackedDataStruct.mat).
%
%   tracks = loadTrackedSeeds()                    all drops, P.seedTracks, 60 fps
%   tracks = loadTrackedSeeds(file, struct('seeds', [5 10], 'fps', 60))
%
% THE DATA. One drop per seed: a single tracked point per video frame (the seed's
% image centroid -- not its CoM, and no attitude). PointsTransformed is the lab
% frame: metres, z vertically UP, seeds released from 10 m. Indices are the video
% frame numbers (60 fps); speed is metres per frame. Each drop starts with 0-3
% frames at rest in the release clamp, so the release is in the record. Some seeds
% glide out of the camera view, so their tracks end early (flagged .partial).
%
% WHAT THIS DOES
%   - time from the frame numbers: t = (frame - releaseFrame)/fps, with the release
%     frame the LAST one at rest, so t = 0 is the moment of release (as in the
%     simulation, which starts from rest);
%   - frames before that are dropped;
%   - zDrop = z(release) - z, the height fallen: the coordinate the comparison
%     with the simulation is aligned on (computePathMetrics windows on it).
% No smoothing here: computePathMetrics smooths simulated and measured paths with
% the same filter.
%
% INPUTS
%   file : (optional) default P.seedTracks from seedPaths.
%   opts : (optional) .seeds (seed numbers; default all), .fps (default 60),
%          .minDrop (m; a track that falls less than this is .partial; default 6).
% OUTPUT
%   tracks : struct array, one per drop: .seedNum, .t (Nx1, s), .pos (Nx3, lab
%            frame, z up, m), .zDrop (Nx1, m), .frames (Nx1), .releaseHeight (m),
%            .fallen (m, total height covered), .partial (logical), .gaps (frames
%            missing inside the track), .fps.

    if nargin < 1 || isempty(file)
        here = fileparts(mfilename('fullpath'));
        while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
        P = seedPaths(here);  file = P.seedTracks;
    end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'fps');     opts.fps     = 60; end
    if ~isfield(opts, 'minDrop'); opts.minDrop = 6;  end
    if ~isfile(file)
        error('loadTrackedSeeds:noFile', 'Tracked-seed file not found: %s', file);
    end

    S = load(file, 'seedTrackedData');   D = S.seedTrackedData;
    nums = [D.seedNum];
    if isfield(opts, 'seeds') && ~isempty(opts.seeds)
        missing = setdiff(opts.seeds, nums);
        if ~isempty(missing)
            error('loadTrackedSeeds:unknownSeed', 'Seed(s) %s not tracked in %s.', mat2str(missing), file);
        end
        idx = arrayfun(@(s) find(nums == s, 1), opts.seeds);
    else
        [~, idx] = sort(nums);
    end

    tracks = struct('seedNum', {}, 't', {}, 'pos', {}, 'zDrop', {}, 'frames', {}, ...
                    'releaseHeight', {}, 'fallen', {}, 'partial', {}, 'gaps', {}, 'fps', {});
    for k = idx(:).'
        d   = D(k);
        pos = d.PointsTransformed;   fr = d.Indices(:);
        % speed(j) is the step from frame j to j+1: the first nonzero one marks the
        % start of motion, so the release frame is j (the last one at rest).
        jMove = find(d.speed > 0, 1);
        if isempty(jMove); jMove = 1; end
        keep = jMove:numel(fr);
        pos  = pos(keep, :);   fr = fr(keep);
        t    = (fr - fr(1)) / opts.fps;
        zDrop = pos(1,3) - pos(:,3);
        fallen = max(zDrop);
        tracks(end+1) = struct('seedNum', d.seedNum, 't', t, 'pos', pos, 'zDrop', zDrop, ...
            'frames', fr, 'releaseHeight', pos(1,3), 'fallen', fallen, ...
            'partial', fallen < opts.minDrop, 'gaps', sum(diff(fr) - 1), 'fps', opts.fps); %#ok<AGROW>
    end
end
