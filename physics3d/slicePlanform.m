function [strips, massParams] = slicePlanform(bsp)
% SLICEPLANFORM  Slice a flat polyshape planform into spanwise strips, and compute
%   the flat-plate mass properties (CoM(t), inertia(t) and its rate) of the wing
%   plus a time-varying point-mass nut.
%
% The geometry half of the shape3d builder (setupSeedShape3D). The builder starts
% from these flat strips and then gives each one its 3D frame (twist, curvature).
%
% PROVENANCE: the computation below is a VERBATIM port of sections 0-6 of the
% frozen planar builder setupSeedShapeAndMass (git tag planar-final), so the
% shape3d seed is bit-identical to what it was when it was built by wrapping that
% builder. Only the input (bsp directly, not seedParamsIn.baseSeedParams) and the
% output packing differ. The planar builder's physics switches and model tag are
% NOT here: they belonged to the planar model.
%
% The input polyshape is defined in a "drawing frame":
%   drawing x -> body z  (spanwise)
%   drawing y -> body x  (chordwise)
% All positions returned use the BODY frame (x=chord, y=normal, z=span).
%
% INPUT
%   bsp : base seed params, with
%     .seedShape       - polyshape of the 2D wing planform (drawing axes above).
%     .nutMass_t       - Nx1 nut mass at each time sample (kg); constant if fixed.
%     .nutPos_t        - 3xN nut position in BODY frame at each time sample (m).
%     .seedDensity     - areal density of the wing (kg/m^2).
%     .seedThickness   - uniform wing thickness (m).
%     .numStrips       - number of spanwise strips.
%     .tSamples        - Nx1 time values for nutMass_t / nutPos_t (s).
%     .plateBoundaries - (OPTIONAL) drawing-x cut positions; default numStrips
%                        equal intervals over the bounding box.
%     .liftMult, .dragMult - (OPTIONAL) per-strip translational lift / drag
%                        multipliers, scalar or 1xM. Default ones.
%
% OUTPUT
%   strips     : .z_body, .chord, .dz, .area, .xgc_body, .zgc_body, .liftMult,
%                .dragMult (each 1xM). Strip geometric centre = the clipped plate's
%                centroid; chord = its bounding-box chord extent; area = its
%                actual clipped area.
%   massParams : .tSamples, .com_t (3xN), .I_G_t (3x3xN), .I_G_dot_t (3x3xN),
%                .M_total. Each plate is a uniform rectangle of its bounding-box
%                size lying flat in the body x-z plane; the nut is a point mass.

% =========================================================================
% 0. UNPACK INPUTS
% =========================================================================

% (bsp is the input: the base seed params, see the header)

seedShape     = bsp.seedShape;       % polyshape in drawing frame
nutMass_t     = bsp.nutMass_t;       % Nx1 nut mass vs time (kg)
nutPos_t      = bsp.nutPos_t;        % 3xN nut position vs time, body frame (m)
aealDensity   = bsp.seedDensity;     % wing areal density (kg/m^2)
thickness     = bsp.seedThickness;   % wing thickness (m)
numStrips     = bsp.numStrips;       % number of aerodynamic strips
tSamples      = bsp.tSamples;        % Nx1 time vector (s)
numT          = numel(tSamples);     % number of time samples

% Optional plate boundaries (drawing-x = body-z direction)
if isfield(bsp, 'plateBoundaries')
    plateBoundaries_drawX = bsp.plateBoundaries;
else
    % Divide the polyshape bounding box evenly into numStrips intervals
    [xv, ~]              = boundingbox(seedShape);   % [xmin xmax] of polyshape
    plateBoundaries_drawX = linspace(xv(1), xv(2), numStrips + 1);
end
numPlates = numel(plateBoundaries_drawX) - 1;   % number of plates / strips

% Wing total area and mass (time-invariant -- only the nut mass varies)
wingArea = area(seedShape);          % m^2
wingMass = aealDensity * wingArea;   % kg


% =========================================================================
% 1. BUILD STRIP GEOMETRY (aerodynamic strips)
%    One strip per plate; the strip width and chord are taken from the plate
%    geometry so strip theory and inertia calculation share the same slicing.
% =========================================================================

% Pre-allocate strip arrays (one entry per plate)
strip_z_body  = zeros(1, numPlates);   % spanwise centre, body z (m)
strip_chord   = zeros(1, numPlates);   % local chord length (m)
strip_dz      = zeros(1, numPlates);   % strip spanwise width (m)
strip_xgc     = zeros(1, numPlates);   % geometric-centre chordwise position, body x (m)
strip_zgc     = zeros(1, numPlates);   % geometric-centre spanwise position,  body z (m)

% Per-plate geometric quantities needed for inertia (computed below)
plate_area        = zeros(1, numPlates);   % plate area (m^2)
plate_centroid_bx = zeros(1, numPlates);   % plate centroid, body x (chordwise, m)
plate_centroid_bz = zeros(1, numPlates);   % plate centroid, body z (spanwise, m)
% Bounding-box dimensions for each plate (used for local inertia approximation)
plate_bbox_dx     = zeros(1, numPlates);   % bounding-box chord extent, body x (m)
plate_bbox_dz     = zeros(1, numPlates);   % bounding-box span  extent, body z (m)

for i = 1 : numPlates

    % --- Spanwise bounds of this plate in drawing-x = body-z ---------------
    z_lo = plateBoundaries_drawX(i);       % left edge  (drawing x = body z)
    z_hi = plateBoundaries_drawX(i + 1);   % right edge

    % --- Clip the polyshape to this spanwise (drawing-x) band ---------------
    % Build a tall rectangle covering the full chordwise extent of the shape
    % at this spanwise position.
    [~, yv] = boundingbox(seedShape);      % [ymin ymax] = chordwise bounds
    clipRect = polyshape( ...
        [z_lo, z_hi, z_hi, z_lo], ...     % drawing-x corners
        [yv(1), yv(1), yv(2), yv(2)] );   % drawing-y corners (full chord)
    platePoly = intersect(seedShape, clipRect);

    if platePoly.NumRegions == 0
        % This slice misses the shape entirely -- skip (area stays zero)
        continue
    end

    % --- Area and centroid of this plate in drawing frame ------------------
    pArea = area(platePoly);               % m^2
    [cx_draw, cy_draw] = centroid(platePoly);
    %   drawing x -> body z,  drawing y -> body x

    % --- Bounding box of the clipped plate polygon -------------------------
    % The bounding box defines the local chord (body x extent) used for both
    % the strip aerodynamics and the local inertia approximation.
    [xv_plate, yv_plate] = boundingbox(platePoly);
    localChord  = yv_plate(2) - yv_plate(1);   % body x extent (chordwise, m)
    localSpan   = xv_plate(2) - xv_plate(1);   % body z extent (spanwise, m)

    % --- Pack into strip arrays (coordinate mapping: draw_x->body_z, draw_y->body_x)
    strip_z_body(i) = (z_lo + z_hi) / 2;   % spanwise centre in body z (m)
    strip_dz(i)     = z_hi - z_lo;          % spanwise width in body z (m)
    strip_chord(i)  = localChord;           % chordwise length in body x (m)

    % Strip GEOMETRIC CENTRE = the plate centroid (both axes). This is a FIXED
    % body-frame aerodynamic reference (mid-chord/centroid), NOT the centre of
    % pressure. The true CoP migrates with angle of attack and is computed each
    % dynamics step in seed6DOFODE (computeAeroCoeffs + computeStripCoP); it is
    % never stored here. The spanwise position is fixed at the strip centroid
    % because strip theory assumes uniform spanwise loading.
    strip_xgc(i) = cy_draw;   % drawing y -> body x (chordwise geometric centre, m)
    strip_zgc(i) = cx_draw;   % drawing x -> body z (spanwise  geometric centre, m)

    % --- Store for inertia calculation below --------------------------------
    plate_area(i)        = pArea;
    plate_centroid_bz(i) = cx_draw;   % drawing x -> body z
    plate_centroid_bx(i) = cy_draw;   % drawing y -> body x
    plate_bbox_dx(i)     = localChord;
    plate_bbox_dz(i)     = localSpan;

end


% =========================================================================
% 2. TIME-VARYING CENTER OF MASS
%    CoM = weighted average of nut position and wing centroid position,
%    computed at each time sample.
%    NOTE: the wing centroid is fixed in the body frame; only the nut mass
%    or position varies in time.
% =========================================================================

% Wing centroid in body frame (time-invariant)
[wingCentroid_drawX, wingCentroid_drawY] = centroid(seedShape);
wingCentroid_body = [wingCentroid_drawY; ...   % drawing y -> body x
                     0;                  ...   % body y = 0 (flat plate)
                     wingCentroid_drawX];      % drawing x -> body z

com_t = zeros(3, numT);   % CoM in body-datum coords at each time sample (m)

for k = 1 : numT
    mNut   = nutMass_t(k);               % nut mass at this time (kg)
    rNut   = nutPos_t(:, k);             % nut position, body frame (3x1, m)
    mTotal = wingMass + mNut;            % total mass at this time (kg)

    % Weighted centroid: (wing contribution + nut contribution) / total mass
    com_t(:, k) = (wingMass * wingCentroid_body + mNut * rNut) / mTotal;
end


% =========================================================================
% 3. TOTAL MASS
%    Assumed constant over time (nut mass fixed). If nutMass_t varies,
%    replace with a vector M_t = wingMass + nutMass_t.
% =========================================================================
M_total = wingMass + nutMass_t(1);   % scalar (kg)


% =========================================================================
% 4. MOMENT OF INERTIA ABOUT THE CoM  (time-varying via parallel-axis theorem)
%
%    I_G(t) = sum_plates [ I_plate_local + m_plate * parallelAxis(d_plate(t)) ]
%             + m_nut(t) * parallelAxis(d_nut(t))
%
%    where d_plate(t) = displacement from plate centroid to CoM at time t.
%
%    I_plate_local is approximated as the inertia of a uniform solid rectangle
%    with the same bounding-box dimensions as the clipped plate polygon.
%    For a rectangle of mass m, chordwise width dx (body x) and spanwise
%    width dz (body z), lying flat in the body x-z plane (body y = 0):
%
%      Ixx_local = (1/12) * m * dz^2          (rotation about body x)
%      Iyy_local = (1/12) * m * (dx^2 + dz^2) (rotation about body y, normal)
%      Izz_local = (1/12) * m * dx^2          (rotation about body z)
%      off-diagonal terms = 0  (rectangle is symmetric about its centroid)
%
%    The off-diagonal terms are zero because the rectangle's centroid coincides
%    with its geometric centre and the axes are aligned with the body frame.
%    This is a good approximation when the plate slice is not strongly skewed.
% =========================================================================

I_G_t = zeros(3, 3, numT);   % 3x3 inertia tensor at each time sample (kg·m^2)

for k = 1 : numT
    I_G_k = zeros(3, 3);    % accumulator for this timestep
    com_k = com_t(:, k);    % CoM at this time (3x1, m)

    % --- Contribution from each plate ---------------------------------------
    for i = 1 : numPlates
        if plate_area(i) == 0
            continue   % plate had no intersection with the shape
        end

        plateMass = aealDensity * plate_area(i);   % plate mass (kg)

        % Plate centroid in body frame (body y = 0, flat plate)
        r_plate_body = [plate_centroid_bx(i); 0; plate_centroid_bz(i)];

        % Displacement from plate centroid to CoM (3x1, m)
        d = r_plate_body - com_k;

        % --- Local inertia of this plate about its own centroid ------------
        % Approximated as a uniform rectangle with the plate bounding-box dims.
        dx = plate_bbox_dx(i);   % chordwise extent, body x (m)
        dz = plate_bbox_dz(i);   % spanwise  extent, body z (m)

        Ixx_local = (1/12) * plateMass * dz^2;           % about body x
        Iyy_local = (1/12) * plateMass * (dx^2 + dz^2);  % about body y (normal)
        Izz_local = (1/12) * plateMass * dx^2;           % about body z

        % Assemble local inertia tensor (no off-diagonal terms for a
        % symmetric rectangle aligned with the body axes)
        I_local = diag([Ixx_local, Iyy_local, Izz_local]);

        % --- Parallel-axis tensor: m * (|d|^2 * I3 - d*d') ----------------
        I_parallel = plateMass * (dot(d, d) * eye(3) - d * d');

        I_G_k = I_G_k + I_local + I_parallel;
    end

    % --- Contribution from the nut (point mass) ------------------------------
    mNut = nutMass_t(k);
    rNut = nutPos_t(:, k);
    d_nut = rNut - com_k;   % displacement from CoM to nut (3x1, m)

    I_G_k = I_G_k + mNut * (dot(d_nut, d_nut) * eye(3) - d_nut * d_nut');

    I_G_t(:, :, k) = I_G_k;
end


% =========================================================================
% 5. NUMERICAL TIME DERIVATIVE OF I_G  (Idot)
%    Central differences in the interior, one-sided at the endpoints.
% =========================================================================

I_G_dot_t = zeros(3, 3, numT);   % d I_G / dt at each time sample (kg·m^2/s)

if numT == 1
    % Single time sample: mass properties are constant, so Idot = 0.
    % I_G_dot_t is already initialised to zeros -- nothing to do.
else
    for k = 1 : numT
        if k == 1
            % Forward difference at the first point
            dt = tSamples(2) - tSamples(1);
            I_G_dot_t(:, :, k) = (I_G_t(:, :, 2) - I_G_t(:, :, 1)) / dt;

        elseif k == numT
            % Backward difference at the last point
            dt = tSamples(end) - tSamples(end - 1);
            I_G_dot_t(:, :, k) = (I_G_t(:, :, end) - I_G_t(:, :, end-1)) / dt;

        else
            % Central difference at interior points
            dt = tSamples(k + 1) - tSamples(k - 1);
            I_G_dot_t(:, :, k) = (I_G_t(:, :, k+1) - I_G_t(:, :, k-1)) / dt;
        end
    end
end


% =========================================================================
% 6. PER-STRIP LIFT / DRAG MULTIPLIERS
%    Static tuning / non-uniform-morphology knobs on the translational lift and
%    drag. Resolve the optional scalar-or-vector inputs into 1xM row vectors
%    (one value per strip); default to all ones (no change) when absent.
% =========================================================================
strip_liftMult = resolveStripMultiplier(bsp, 'liftMult', numPlates);
strip_dragMult = resolveStripMultiplier(bsp, 'dragMult', numPlates);

% =========================================================================
% 7. PACK THE OUTPUTS
% =========================================================================
strips.z_body   = strip_z_body;   % 1xM, body-z centres (m)
strips.chord    = strip_chord;    % 1xM, chord lengths (m)
strips.dz       = strip_dz;       % 1xM, spanwise widths (m)
strips.area     = plate_area;     % 1xM, ACTUAL clipped plate area (m^2), so per-strip
                                  % masses conserve the wing mass on tapered planforms
strips.xgc_body = strip_xgc;      % 1xM, geometric-centre chordwise body x (m)
strips.zgc_body = strip_zgc;      % 1xM, geometric-centre spanwise  body z (m)
strips.liftMult = strip_liftMult; % 1xM, translational-lift multipliers
strips.dragMult = strip_dragMult; % 1xM, drag multipliers

massParams.tSamples  = tSamples;  % Nx1 (s)
massParams.com_t     = com_t;     % 3xN, body-datum (m)
massParams.I_G_t     = I_G_t;     % 3x3xN (kg m^2)
massParams.I_G_dot_t = I_G_dot_t; % 3x3xN (kg m^2/s)
massParams.M_total   = M_total;   % scalar (kg)

end

% =========================================================================
% LOCAL: resolve an optional scalar-or-vector strip multiplier to 1xM
% =========================================================================
function multVec = resolveStripMultiplier(bsp, fieldName, numStrips)
% RESOLVESTRIPMULTIPLIER  Turn an optional lift/drag multiplier input into a
%   1xM row vector, one value per strip.
%
% INPUTS
%   bsp       : base seed-params struct (may or may not contain fieldName).
%   fieldName : multiplier field name, e.g. 'liftMult' or 'dragMult'.
%   numStrips : number of strips M the output must span.
%
% OUTPUT
%   multVec   : 1xM row vector. All ones if the field is absent/empty; a scalar
%               is broadcast to every strip; a length-M vector is used as-is
%               (reshaped to a row). Any other length is an error.

    if ~isfield(bsp, fieldName) || isempty(bsp.(fieldName))
        multVec = ones(1, numStrips);       % default: no scaling
        return
    end

    raw = bsp.(fieldName)(:).';              % force a row vector

    if isscalar(raw)
        multVec = raw * ones(1, numStrips);  % broadcast one value to all strips
    elseif numel(raw) == numStrips
        multVec = raw;                       % explicit per-strip profile
    else
        error('slicePlanform:badMultiplierLength', ...
              '%s must be a scalar or a length-%d vector (got length %d).', ...
              fieldName, numStrips, numel(raw));
    end
end
