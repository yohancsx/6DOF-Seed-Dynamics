function visualizeSeedLocalVels(seedParams, seedLocalVels, rot, pos, opts)
% VISUALIZESEEDLOCALVELS  Draw seed local velocity vectors at a given pose.
%
% Plots velocity vectors (from computeSeedLocalVel) as scaled lines with a dot
% at the tip. Can draw the per-strip velocities at each geometric centre, or
% a single averaged vector at the CoM. Uses the same Y-up plot convention as
% visualizeSeedShape and does NOT clear the figure, so it layers on top of a
% seed-shape plot.
%
% INPUTS
%   seedParams    : struct from setupSeedShapeAndMass (strip geometry, CoM).
%   seedLocalVels : struct from computeSeedLocalVel with 3xM body-frame fields
%                   .totalVel, .windVel, .projectedWind.
%   rot           : 3x3 rotation matrix, body->world. Pass eye(3) for none.
%   pos           : 3x1 world-frame position of the body-datum origin (m).
%   opts          : (optional) struct of toggles, colours, and scale. Defaults:
%     .tIndex         - time index for the CoM lookup       (default: 1)
%     .mode           - 'geoCenter' (per strip) or 'com' (average) (default: 'geoCenter')
%     .showTotal      - draw totalVel vectors               (default: true)
%     .showWind       - draw windVel vectors                (default: false)
%     .showProjected  - draw projectedWind vectors          (default: false)
%     .colorTotal     - totalVel colour       (default: [0.20 0.40 0.80] blue)
%     .colorWind      - windVel colour        (default: [0.85 0.40 0.10] orange)
%     .colorProjected - projectedWind colour  (default: [0.55 0.15 0.60] purple)
%     .scale          - arrow length per unit speed, m/(m/s)  (default: 0.03)
%     .dotSize        - marker size of the tip dot            (default: 6)
%     .lineWidth      - arrow line width                      (default: 1.5)

% -------------------------------------------------------------------------
% 0. DEFAULT OPTIONS
% -------------------------------------------------------------------------
if nargin < 5 || isempty(opts)
    opts = struct();
end

opts = setDefault(opts, 'tIndex',        1);
opts = setDefault(opts, 'mode',          'geoCenter');
opts = setDefault(opts, 'showTotal',     true);
opts = setDefault(opts, 'showWind',      false);
opts = setDefault(opts, 'showProjected', false);
opts = setDefault(opts, 'colorTotal',    [0.20 0.40 0.80]);   % blue
opts = setDefault(opts, 'colorWind',     [0.85 0.40 0.10]);   % orange
opts = setDefault(opts, 'colorProjected',[0.55 0.15 0.60]);   % purple
opts = setDefault(opts, 'scale',         0.03);
opts = setDefault(opts, 'dotSize',       6);
opts = setDefault(opts, 'lineWidth',     1.5);

% -------------------------------------------------------------------------
% 1. WORLD->PLOT MAPPING  (Y-up, matches visualizeSeedShape)
%    plot X = world X,  plot Y = world Z,  plot Z = world Y (vertical)
% -------------------------------------------------------------------------
toPlot = @(w) [w(1,:); w(3,:); w(2,:)];

% -------------------------------------------------------------------------
% 2. ARROW BASE POSITIONS (body frame)
% -------------------------------------------------------------------------
xgc = seedParams.strips.xgc_body;   % 1xM chordwise geometric centre, body x (m)
zgc = seedParams.strips.zgc_body;   % 1xM spanwise  geometric centre, body z (m)
numStrips = numel(xgc);
geoCenterBody = [xgc; zeros(1, numStrips); zgc];   % 3xM, body y = 0 (flat plate)

% CoM in body-datum coords (used as the base point in 'com' mode)
comBody = seedParams.massParams.com_t(:, opts.tIndex);   % 3x1

% -------------------------------------------------------------------------
% 3. BUILD LIST OF VELOCITY SETS TO DRAW  {data(3xM), colour}
% -------------------------------------------------------------------------
sets = {};
if opts.showTotal
    sets(end+1, :) = {seedLocalVels.totalVel,      opts.colorTotal};
end
if opts.showWind
    sets(end+1, :) = {seedLocalVels.windVel,       opts.colorWind};
end
if opts.showProjected
    sets(end+1, :) = {seedLocalVels.projectedWind, opts.colorProjected};
end

% -------------------------------------------------------------------------
% 4. DRAW
% -------------------------------------------------------------------------
for s = 1 : size(sets, 1)
    velAll = sets{s, 1};   % 3xM, body frame
    col    = sets{s, 2};

    if strcmpi(opts.mode, 'com')
        % Single averaged vector at the CoM
        velAvg = mean(velAll, 2);   % 3x1
        drawArrow(comBody, velAvg, rot, pos, toPlot, opts, col);
    else
        % One vector per strip at its geometric centre
        for i = 1 : numStrips
            drawArrow(geoCenterBody(:, i), velAll(:, i), rot, pos, toPlot, opts, col);
        end
    end
end

end   % visualizeSeedLocalVels


% =========================================================================
% LOCAL HELPER: draw one scaled velocity arrow (line + tip dot)
% =========================================================================
function drawArrow(baseBody, velBody, rot, pos, toPlot, opts, col)
% Base is a POSITION (rotated and translated); velocity is a free VECTOR
% (rotated only, never translated). Both are then mapped to plot axes.
    baseWorld = rot * baseBody + pos;                 % 3x1 position
    tipWorld  = baseWorld + opts.scale * (rot * velBody);   % 3x1 position

    bp = toPlot(baseWorld);
    tp = toPlot(tipWorld);

    % Arrow shaft
    plot3([bp(1), tp(1)], [bp(2), tp(2)], [bp(3), tp(3)], ...
          '-', 'Color', col, 'LineWidth', opts.lineWidth);
    hold on;

    % Tip dot (sits at the base when the vector is zero, e.g. projected wind
    % for purely spanwise motion -- a useful visual cue)
    plot3(tp(1), tp(2), tp(3), ...
          'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', opts.dotSize);
    hold on;
end


% =========================================================================
% LOCAL HELPER: set default struct field if absent
% =========================================================================
function s = setDefault(s, field, value)
    if ~isfield(s, field)
        s.(field) = value;
    end
end




% function visualizeSeedLocalVels(seedParams,seedLocalVels,rot, pos, opts)
%     %visualize the seed local velocities
%     %the options should contain the possible colors for each velocity, and
%     %which velocities to plot. This code can either plot the average of the
%     %seedLocalVels at the COM or the seed local vels at each center of
%     %pressure. Then, there should be a choice of which velocity to plot
%     %based on opts. Plot the velocity as a line with a small dot at the end
%     %of a given input size in opts (the line and dot should have a scale). 
% end