function seedParamsFull = setupSeedShapeAndMass(seedParamsIn)
% SETUPSEEDSHAPEANDMASS  Populate all geometry and mass fields of seedParams
%   from a 2D polyshape definition of the seed wing + a time-varying nut mass.
%
% The input polyshape is defined in a "drawing frame":
%   drawing x -> body z  (spanwise)
%   drawing y -> body x  (chordwise)
% All positions stored in seedParamsFull use the BODY frame (x=chord, y=normal,
% z=span), consistent with the dynamics skeleton.
%
% INPUT
%   seedParamsIn : struct with sub-struct seedParamsIn.baseSeedParams containing:
%     .seedShape       - polyshape object defining the 2D wing planform.
%                        Drawing axes: x = spanwise (-> body z),
%                                      y = chordwise (-> body x).
%     .nutMass_t       - Nx1 vector of nut masses at each time sample (kg).
%                        Use a constant vector if the nut mass is fixed.
%     .nutPos_t        - 3xN array of nut positions in BODY frame at each
%                        time sample (m). Each column is [x;y;z] at that time.
%     .seedDensity     - areal density of the seed wing (kg/m^2).
%     .seedThickness   - uniform wing thickness (m), used for volume / buoyancy.
%     .numStrips       - number of spanwise strips to divide the shape into.
%     .tSamples        - Nx1 vector of time values corresponding to nutMass_t
%                        and nutPos_t (s).
%     .plateBoundaries - (OPTIONAL) vector of drawing-x values at which to cut
%                        the shape into discrete plates for inertia calculation.
%                        If absent, the shape is divided into numStrips equal
%                        spanwise intervals.
%
% OUTPUT
%   seedParamsFull : copy of seedParamsIn with the following fields added:
%     .strips.z_body   - 1xM vector of strip spanwise centres in body frame (m).
%     .strips.chord    - 1xM vector of local chord lengths at each strip (m).
%     .strips.dz       - 1xM vector of strip spanwise widths (m).
%     .strips.xgc_body - 1xM vector of strip GEOMETRIC-CENTRE x positions
%                        in body frame (chordwise), = the strip centroid x.
%                        This is a fixed aerodynamic reference (mid-chord), NOT
%                        the centre of pressure.
%     .strips.zgc_body - 1xM vector of strip GEOMETRIC-CENTRE z positions
%                        in body frame (spanwise), = the strip centroid z.
%     .massParams.tSamples - copy of tSamples (s).
%     .massParams.com_t    - 3xN array, CoM position in body-DATUM coords at
%                            each time sample (m).
%     .massParams.I_G_t    - 3x3xN array, inertia tensor about CoM at each
%                            time sample (kg·m^2).
%     .massParams.I_G_dot_t- 3x3xN array, numerical time-derivative of I_G (kg·m^2/s).
%     .massParams.M_total  - total mass of the seed (nut + wing), scalar (kg).
%                            Assumed constant; if the nut mass varies, this
%                            becomes a vector -- extend as needed.
%     .baseSeedParams      - pass-through of the original input sub-struct.

% =========================================================================
% 0. UNPACK INPUTS
% =========================================================================

bsp = seedParamsIn.baseSeedParams;   % short alias for readability

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
% 6. PACK ALL OUTPUTS
% =========================================================================

seedParamsFull = seedParamsIn;   % carry through all original inputs

% --- Strip geometry (used by the aerodynamic strip loop) -----------------
seedParamsFull.strips.z_body   = strip_z_body;   % 1xM, body-z centres (m)
seedParamsFull.strips.chord    = strip_chord;    % 1xM, chord lengths (m)
seedParamsFull.strips.dz       = strip_dz;       % 1xM, spanwise widths (m)
seedParamsFull.strips.xgc_body = strip_xgc;      % 1xM, geometric-centre chordwise body x (m)
seedParamsFull.strips.zgc_body = strip_zgc;      % 1xM, geometric-centre spanwise  body z (m)

% --- Time-varying mass properties (queried by getMassProperties) ---------
seedParamsFull.massParams.tSamples   = tSamples;     % Nx1 (s)
seedParamsFull.massParams.com_t      = com_t;         % 3xN, body-datum (m)
seedParamsFull.massParams.I_G_t      = I_G_t;         % 3x3xN (kg·m^2)
seedParamsFull.massParams.I_G_dot_t  = I_G_dot_t;     % 3x3xN (kg·m^2/s)
seedParamsFull.massParams.M_total    = M_total;        % scalar (kg)

% --- Pass through original sub-struct ------------------------------------
seedParamsFull.baseSeedParams = bsp;

end

%function seedParamsFull = setupSeedShapeAndMass(seedParamsIn)
    %take the input seed params (seedParamsIn) and fully populate the
    %seedParams. Made for a 2D input seed. All values input in kg, m etc.
    %INPUT: seedParamsIn contains inside struct variable
    %seedParamsIn.baseSeedParams
    % - mass and location of the seed nut over time (idealized as a point mass)
    % - a 2D area defined by a polyshape which defines the seed shape, the
    % horizontal (x direction) of this shape is the body frame z, and the
    % vertical direction is the body frame x
    % - seed density, which is the density of the seed 2D area
    % - seed thickness (assumed to be constant)
    % - number of seed plates
    % - (optional) plateBoundaries - a vector of x values where to cut the
    % seed into plates

    %OUTPUT: should populate the seedParams struct with the necessary data
    %to run the simulation
    % - the COM position over time calculated by using the mass and
    % location of the seed nut as well as the polyshape
    % - where the seed strips are (strips for the thin plate theory), break
    % the polyshape up into strips uniformly spaced along the input x axis
    % (this is really the z axis in the body frame, so make sure to set it
    % up as such with the input y axis being the new x)
    % - the centers of pressures of the strips
    % - the total moment of inertia of the seed over time at the center of mass,
    % and it's derivative
    % - all the inputs in seedParamsIn

    %parse the inputs

    %if there is a variable in the seedParamsIn called plateBoundaries,
    %then use those as the boundaries for the thin plates, otherwise, set
    %the plate boundaries by dividing the input polyshape evenly along the
    %x axis 

    %for each plate, determine the center of area, and set the center of
    %pressure to that point (X and Y position)

    %for each plate, determine the individual moment of inertia of each
    %plate. Do this by setting a rectangular bounding box for each plate
    %and computing the moment of inertia

    %use the paralell axis theorem to sum the moments of inertia about the
    %COM (do this for each timestep). make sure to account for the nut as a
    %point mass. Since all the plate Inertia axes should be paralell, this
    %should be fairly simple.

    %take the derivative numerically to determine Idot

    %store all values and return
%end