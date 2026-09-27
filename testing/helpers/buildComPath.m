function [tDense, nutPos, eventT, dwellInt] = buildComPath(posList, dwellList, moveTime, dt)
% BUILDCOMPATH  Nut-position waypoints (+ dwells) -> dense, pchip-smoothed path.
%
% Builds a strictly-increasing waypoint schedule -- hold each posList column for
% its dwell, then move to the next over moveTime -- and interpolates it onto a
% dense grid with pchip (shape-preserving: flat dwells stay flat, moves ease in
% and out, no overshoot past the commanded limits).
%
% INPUTS
%   posList   : 3xK nut-position waypoints, body [x;y;z] (m).
%   dwellList : 1xK time held at each waypoint (s). A zero dwell is allowed.
%   moveTime  : transition time between consecutive waypoints (s).
%   dt        : sample step of the returned dense path (s).
%
% OUTPUTS
%   tDense   : 1xN time grid, 0 .. total schedule duration (s).
%   nutPos   : 3xN interpolated nut position (m).
%   eventT   : 1xK arrival time of each waypoint (= start of its dwell), for
%              marking the CoM-move instants on downstream plots/animations.
%   dwellInt : Kx2 [start end] of each dwell, for classifying the SETTLED
%              behaviour at a waypoint rather than the transient into it.
%
% This is the shared version of a routine that previously existed only as a
% local function inside testing/planar/runComMovementTest.m and
% testing/baselines/generateModelBaseline.m; those copies are untouched (the
% baseline generator is pinned by its snapshots), but new callers use this one.

    K = size(posList, 2);
    wpT = 0;                  % waypoint times
    wpP = posList(:, 1);      % waypoint positions (3 x .)
    tcur = 0;
    eventT   = zeros(1, K);   % arrival time of each waypoint (dwell start)
    dwellInt = zeros(K, 2);   % [start end] of each dwell
    for i = 1:K
        eventT(i)     = tcur;
        dwellInt(i,1) = tcur;
        if dwellList(i) > 0                       % hold at waypoint i
            tcur = tcur + dwellList(i);
            wpT(end+1)    = tcur;             %#ok<AGROW>
            wpP(:, end+1) = posList(:, i);    %#ok<AGROW>
        end
        dwellInt(i,2) = tcur;
        if i < K                                  % move to waypoint i+1
            tcur = tcur + moveTime;
            wpT(end+1)    = tcur;             %#ok<AGROW>
            wpP(:, end+1) = posList(:, i+1);  %#ok<AGROW>
        end
    end

    tDense = 0:dt:wpT(end);
    if tDense(end) < wpT(end); tDense(end+1) = wpT(end); end   % include the endpoint

    nutPos = zeros(3, numel(tDense));
    for a = 1:3
        nutPos(a, :) = interp1(wpT, wpP(a, :), tDense, 'pchip');
    end
end
