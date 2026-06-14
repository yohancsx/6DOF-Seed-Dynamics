function visualizeSeedShape(seedParams, rot, pos, opts)
% VISUALIZESEEDSHAPE  Plot the seed geometry in 3D at a given pose.
%
% Renders the seed wing as a filled surface, with optional overlays for
% strip plate boundaries, per-plate CoM markers, per-strip CoP markers,
% the seed CoM, and the nut. Does NOT clear the figure or call hold off,
% so it can be layered into any existing axes.
%
% INPUTS
%   seedParams : struct output of setupSeedShapeAndMass.
%   rot        : 3x3 rotation matrix, body->world (columns are body axes in
%                world coords). Pass eye(3) for no rotation.
%   pos        : 3x1 world-frame position of the body-frame datum origin (m).
%   opts       : (optional) struct of display toggles and colours. Any
%                omitted field falls back to the default listed below.
%     .tIndex         - index into time-varying arrays (default: 1)
%     .showPlate      - draw the wing surface patch      (default: true)
%     .showPlateEdges - draw strip boundary lines        (default: true)
%     .showPlateCom   - draw per-plate centroid markers  (default: false)
%     .showCop        - draw per-strip CoP markers       (default: false)
%     .showCom        - draw total seed CoM marker       (default: true)
%     .showNut        - draw the nut position marker     (default: true)
%     .colorPlate     - wing surface face color          (default: [0.18 0.55 0.18])
%     .colorEdge      - strip boundary line color        (default: [0.10 0.35 0.10])
%     .colorPlateCom  - per-plate CoM marker color       (default: [0.80 0.20 0.20])
%     .colorCom       - total CoM marker color           (default: [0.85 0.10 0.10])
%     .colorCop       - CoP marker color                 (default: [0.20 0.40 0.80])
%     .colorNut       - nut marker color                 (default: [0.10 0.10 0.10])
%     .plateAlpha     - wing surface transparency 0-1    (default: 0.7)
%     .markerSize     - size of all point markers        (default: 10)

% -------------------------------------------------------------------------
% 0. DEFAULT OPTIONS
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opts)
    opts = struct();
end

% Toggles
opts = setDefault(opts, 'tIndex',         1);
opts = setDefault(opts, 'showPlate',      true);
opts = setDefault(opts, 'showPlateEdges', true);
opts = setDefault(opts, 'showPlateCom',   false);
opts = setDefault(opts, 'showCop',        false);
opts = setDefault(opts, 'showCom',        true);
opts = setDefault(opts, 'showNut',        true);

% Colours
opts = setDefault(opts, 'colorPlate',    [0.18 0.55 0.18]);   % mid green
opts = setDefault(opts, 'colorEdge',     [0.10 0.35 0.10]);   % dark green
opts = setDefault(opts, 'colorPlateCom', [0.80 0.20 0.20]);   % red
opts = setDefault(opts, 'colorCom',      [0.85 0.10 0.10]);   % bright red
opts = setDefault(opts, 'colorCop',      [0.20 0.40 0.80]);   % blue
opts = setDefault(opts, 'colorNut',      [0.10 0.10 0.10]);   % black

% Style
opts = setDefault(opts, 'plateAlpha', 0.7);
opts = setDefault(opts, 'markerSize', 10);

% -------------------------------------------------------------------------
% 1. UNPACK seedParams FIELDS
% -------------------------------------------------------------------------
bsp      = seedParams.baseSeedParams;
strips   = seedParams.strips;
mp       = seedParams.massParams;
tIdx     = opts.tIndex;

% Wing polyshape in drawing frame (draw_x=span=body_z, draw_y=chord=body_x)
wingPoly  = bsp.seedShape;

% CoM at this time step, body-datum coords (3x1)
com_body  = mp.com_t(:, tIdx);

% Nut position and mass at this time step
nutPos_body = bsp.nutPos_t(:, tIdx);   % body frame (3x1, m)
nutMass     = bsp.nutMass_t(tIdx);     % kg

% Strip geometry
z_body  = strips.z_body;    % 1xM spanwise centres, body z (m)
chord   = strips.chord;     % 1xM chord lengths,    body x (m)
dz      = strips.dz;        % 1xM spanwise widths,  body z (m)
xcp     = strips.xcp_body;  % 1xM CoP chordwise,    body x (m)
zcp     = strips.zcp_body;  % 1xM CoP spanwise,     body z (m)
numStrips = numel(z_body);

% -------------------------------------------------------------------------
% 2. HELPERS
%    toWorld  : body frame -> world frame  (physics convention)
%    toPlot   : world frame -> plot axes   (display convention)
%
%    Display convention: body/world Y is plate-normal, which aligns with
%    gravity (down = -Y). To show this on the vertical axis we remap:
%        plot X = world X  (chordwise when rot = eye)
%        plot Y = world Z  (spanwise  when rot = eye)
%        plot Z = world Y  (normal, vertical -- up is +Y)
%    All plot3/patch calls use toPlot(toWorld(bodyPts)) so the mapping is
%    applied consistently in one place.
% -------------------------------------------------------------------------
toWorld = @(bodyPts) rot * bodyPts + pos;            % body  -> world, 3xN
toPlot  = @(w) [w(1,:); w(3,:); w(2,:)];            % world -> plot axes, 3xN

% -------------------------------------------------------------------------
% 3. WING SURFACE PATCH
%    The polyshape lives in drawing frame (draw_x = body_z, draw_y = body_x,
%    body_y = 0). Extrude a zero-thickness flat patch in body_y = 0.
% -------------------------------------------------------------------------
if opts.showPlate
    % Extract polyshape boundary vertices (drawing frame)
    [vx_draw, vy_draw] = boundary(wingPoly);   % spanwise, chordwise (drawing)

    % Map to body frame: body_x = draw_y, body_y = 0, body_z = draw_x
    bodyVerts = [vy_draw'; ...                  % body x (chordwise)
                 zeros(1, numel(vx_draw));  ...  % body y = 0 (flat plate)
                 vx_draw'];                      % body z (spanwise)

    % Rotate and translate to world frame, then remap to plot axes
    p = toPlot(toWorld(bodyVerts));   % 3 x nVerts, plot convention

    patch(p(1,:), p(2,:), p(3,:), ...
          opts.colorPlate, ...
          'FaceAlpha', opts.plateAlpha, ...
          'EdgeColor', 'none');
    hold on;
end

% -------------------------------------------------------------------------
% 4. STRIP PLATE BOUNDARY EDGES
%    Draw a chordwise line at each strip boundary (left edge of each strip,
%    plus the far right edge of the last strip).
% -------------------------------------------------------------------------
if opts.showPlateEdges
    % Collect all unique boundary z values (body z)
    zEdges = [z_body - dz/2,  z_body(end) + dz(end)/2];   % M+1 edges

    % For each edge, find the chord extent at that z by interpolating the
    % chord array. For a rectangular seed all chords are equal; for a
    % tapered one this gives the correct local chord.
    % Use nearest-strip chord as the edge chord (adequate for visualisation).
    for e = 1 : numel(zEdges)
        ze = zEdges(e);
        % Nearest strip index
        [~, iNear] = min(abs(z_body - ze));
        localChord = chord(iNear);

        % Chordwise endpoints at this spanwise edge (body frame)
        x_lo = -localChord / 2;
        x_hi =  localChord / 2;

        edgePts_body = [x_lo, x_hi; ...   % body x
                        0,    0;    ...   % body y
                        ze,   ze];        % body z

        p = toPlot(toWorld(edgePts_body));
        plot3(p(1,:), p(2,:), p(3,:), ...
              '-', 'Color', opts.colorEdge, 'LineWidth', 1);
        hold on;
    end
end

% -------------------------------------------------------------------------
% 5. PER-PLATE CENTROID MARKERS
% -------------------------------------------------------------------------
if opts.showPlateCom
    for i = 1 : numStrips
        % Strip centroid: chordwise centre = 0 (symmetric), spanwise = z_body
        centroid_body = [0; 0; z_body(i)];
        p = toPlot(toWorld(centroid_body));
        plot3(p(1), p(2), p(3), ...
              'o', 'Color', opts.colorPlateCom, ...
              'MarkerFaceColor', opts.colorPlateCom, ...
              'MarkerSize', opts.markerSize * 0.7);
        hold on;
    end
end

% -------------------------------------------------------------------------
% 6. PER-STRIP CENTRE-OF-PRESSURE MARKERS
% -------------------------------------------------------------------------
if opts.showCop
    for i = 1 : numStrips
        cop_body  = [xcp(i); 0; zcp(i)];
        p = toPlot(toWorld(cop_body));
        plot3(p(1), p(2), p(3), ...
              'd', 'Color', opts.colorCop, ...
              'MarkerFaceColor', opts.colorCop, ...
              'MarkerSize', opts.markerSize * 0.7);
        hold on;
    end
end

% -------------------------------------------------------------------------
% 7. TOTAL SEED CENTER OF MASS
% -------------------------------------------------------------------------
if opts.showCom
    com_world = toPlot(toWorld(com_body));
    plot3(com_world(1), com_world(2), com_world(3), ...
          'o', 'Color', opts.colorCom, ...
          'MarkerFaceColor', opts.colorCom, ...
          'MarkerSize', opts.markerSize);
    hold on;
end

% -------------------------------------------------------------------------
% 8. NUT MARKER
% -------------------------------------------------------------------------
if opts.showNut && nutMass > 0
    nut_world = toPlot(toWorld(nutPos_body));
    plot3(nut_world(1), nut_world(2), nut_world(3), ...
          's', 'Color', opts.colorNut, ...
          'MarkerFaceColor', opts.colorNut, ...
          'MarkerSize', opts.markerSize);
    hold on;
end

end   % visualizeSeedShape


% =========================================================================
% LOCAL HELPER
% =========================================================================
function s = setDefault(s, field, value)
% SETDEFAULT  Set s.field = value if the field is absent from struct s.
    if ~isfield(s, field)
        s.(field) = value;
    end
end

% function visualizeSeedShape(seedParams,rot,pos, opts)
%     %function to visualize the seed shape. Has the rotation and position of
%     %the seed as well. Optionally has the plotting params like if the nut
%     %should be visualized, if the COM of each seed plate should be
%     %visualized, if each plate edges should be visualized, if the COPs of
%     %the seed plates should be realized, etc. The colors can also be
%     %optionally set as well. In this case, the seed plotting may be used in
%     %tandem with other vidualizations, so don't clear the figure and make
%     %it easy to incorporate with other plotting code.
% 
%     %plot the seed plate, seed should be green, nut should be black, COM
%     %should be red and other colors can be decided.
% end