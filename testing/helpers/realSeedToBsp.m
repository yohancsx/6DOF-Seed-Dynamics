function [bsp, meta] = realSeedToBsp(o, parts, opts)
% REALSEEDTOBSP  One measured seed (seedMassModel.mat) -> base seed params with an
%   explicit strip geometry and mass model (bsp.stripGeom), for buildSeedParams.
%
% Called by loadRealSeeds for each seed; separate so testRealSeedLoader can feed it
% a synthetic flat plate and compare with the polyshape builder.
%
% THE DATA (CT-fitted, mm / g / g mm^2, MESH axes: x = span, y = chord, z = normal)
%   o.plates.plates(i)   12 flat trapezoids: frame e1 (span) / e2 (chord) / n
%                        (normal), right-handed (e2 x n = e1); corners (4x3); length
%                        (along e1); area; twistDeg, dihedralDeg; thickness.
%   parts                the part table rows for this seed (TP in the .mat): one row
%                        per plate (body "plate") and one for the nut (body "nut"),
%                        each with mass_g (CALIBRATED), predMass_g (from the
%                        assumed densities), centroidLumped_*, and the inertia
%                        I*_lumped_gmm2 about the part's own centroid. Summed with
%                        parallel-axis terms they reproduce o.calibrated (mass, CoM,
%                        inertia) exactly.
%   o.measuredMass, o.calibration (= measured / predicted mass).
%
% WHAT GOES WHERE IN THE MODEL
%   strips      the aerodynamics. One strip per plate (or opts.subStrips per plate,
%               cut along e1, which is exact for a flat trapezoid): centre at the
%               strip's AREA centroid, dz = its length along e1, chord = area/dz
%               (so sum(chord.*dz) is the plate area), frame from the plate.
%   massParams  the rigid body: total mass, CoM and inertia about the CoM, from the
%               parts. Constant in time. The strips carry NO mass of their own here.
%   Both are in BODY axes (x = chord, y = normal, z = span), from the mesh by the
%   cyclic permutation P = [0 1 0; 0 0 1; 1 0 0] (a proper rotation):
%       positions and frame vectors  v_body = P v_mesh
%       inertia tensors              I_body = P I_mesh P'
%   The plates and the parts share the mesh origin, so the CoM lands in the same
%   datum as the strip centres -- which the added mass and its off-centre coupling
%   rely on (getAddedMass measures strip lever arms from the CoM).
%
% CALIBRATION (opts.calibration)
%   'both'    nut and wing densities scaled together to the measured mass: the
%             parts' mass_g as given (= o.calibrated).
%   'nutOnly' wing at its predicted (assumed-density) mass; the nut takes the rest,
%             m_nut = measured - sum(wing predMass). Each part keeps its centroid,
%             and its own inertia scales with its mass (same distribution).
%
% INPUTS
%   o     : one element of the .mat's `out` struct array.
%   parts : table rows of the .mat's `TP` for this seed.
%   opts  : .calibration ('both' | 'nutOnly'), .subStrips (positive integer).
% OUTPUTS
%   bsp   : base seed params: .stripGeom (.strips, .massParams, .source) plus the
%           fields other code reads -- .seedShape (planform polyshape, drawX = body
%           z, drawY = body x, for the renderers), .seedThickness (area-weighted
%           plate thickness; edge-drag area), .seedDensity (wing areal density),
%           .numStrips, .tSamples = 0, .nutPos_t / .nutMass_t (the nut, for the
%           renderers; NOT used for mass, which comes from stripGeom).
%   meta  : per-seed numbers for reports (masses, span, area, AR, loading, ...).

    P   = [0 1 0; 0 0 1; 1 0 0];          % mesh -> body
    mm  = 1e-3;   g = 1e-3;   gmm2 = 1e-9;
    nS  = opts.subStrips;
    pl  = o.plates.plates;   nP = numel(pl);

    % =====================================================================
    % 1. STRIPS from the plates
    % =====================================================================
    M = nP * nS;
    chord = zeros(1,M);  dz = chord;  area = chord;  ctr = zeros(3,M);
    cDir = zeros(3,M);   nDir = cDir;  sDir = cDir;  twist = chord;  dihed = chord;
    thick = chord;
    k = 0;
    for i = 1:nP
        p  = pl(i);
        e1 = p.e1(:);  e2 = p.e2(:);  n = p.n(:);
        C  = p.corners;                                  % 4x3, mesh
        % inboard / outboard edges by position along e1, each edge's ends by e2
        [~, o1] = sort(C * e1);   inE = C(o1(1:2), :);   outE = C(o1(3:4), :);
        [~, a] = sort(inE * e2);  inE  = inE(a, :);
        [~, b] = sort(outE * e2); outE = outE(b, :);
        for j = 1:nS
            t0 = (j-1)/nS;   t1 = j/nS;
            A0 = inE(1,:) + t0*(outE(1,:) - inE(1,:));   A1 = inE(1,:) + t1*(outE(1,:) - inE(1,:));
            B0 = inE(2,:) + t0*(outE(2,:) - inE(2,:));   B1 = inE(2,:) + t1*(outE(2,:) - inE(2,:));
            [aT1, cT1] = tri(A0, B0, B1);   [aT2, cT2] = tri(A0, B1, A1);
            k = k + 1;
            area(k)  = (aT1 + aT2) * mm^2;
            ctr(:,k) = P * ((aT1*cT1 + aT2*cT2) / (aT1 + aT2)).' * mm;
            dz(k)    = p.length / nS * mm;
            chord(k) = area(k) / dz(k);
            cDir(:,k) = P*e2;   nDir(:,k) = P*n;   sDir(:,k) = P*e1;
            twist(k) = deg2rad(p.twistDeg);   dihed(k) = deg2rad(p.dihedralDeg);
            thick(k) = p.thickness * mm;
        end
    end
    s.z_body    = ctr(3,:);
    s.chord     = chord;
    s.dz        = dz;
    s.area      = area;
    s.xgc_body  = ctr(1,:);
    s.zgc_body  = ctr(3,:);
    s.liftMult  = ones(1,M);
    s.dragMult  = ones(1,M);
    s.ygc_body  = ctr(2,:);
    s.chordDir  = cDir;
    s.normalDir = nDir;
    s.spanDir   = sDir;
    s.twist     = twist;
    s.dihedral  = dihed;

    % =====================================================================
    % 2. MASS from the parts (about the CoM), with the chosen calibration
    % =====================================================================
    isNut = parts.body == "nut";
    mCal  = parts.mass_g;
    mPred = parts.predMass_g;
    switch opts.calibration
        case 'both'
            mUse = mCal;
        case 'nutOnly'
            mUse = mPred;
            mUse(isNut) = o.measuredMass - sum(mPred(~isNut));
            if any(mUse(isNut) <= 0)
                error('realSeedToBsp:nutOnlyNegative', ...
                      '%s: the predicted wing mass exceeds the measured mass; nutOnly impossible.', o.name);
            end
        otherwise
            error('realSeedToBsp:badCalibration', 'opts.calibration must be ''both'' or ''nutOnly''.');
    end
    cen = parts{:, {'centroidLumped_x','centroidLumped_y','centroidLumped_z'}};
    Mtot = sum(mUse);
    com  = (mUse.' * cen) / Mtot;
    I    = zeros(3);
    for r = 1:height(parts)
        Ii = [parts.Ixx_lumped_gmm2(r) parts.Ixy_lumped_gmm2(r) parts.Ixz_lumped_gmm2(r);
              parts.Ixy_lumped_gmm2(r) parts.Iyy_lumped_gmm2(r) parts.Iyz_lumped_gmm2(r);
              parts.Ixz_lumped_gmm2(r) parts.Iyz_lumped_gmm2(r) parts.Izz_lumped_gmm2(r)];
        Ii = Ii * (mUse(r) / mCal(r));                  % own inertia scales with mass
        d  = cen(r,:) - com;
        I  = I + Ii + mUse(r) * (dot(d,d)*eye(3) - d.'*d);
    end
    mpar.tSamples  = 0;
    mpar.com_t     = P * com.' * mm;
    mpar.I_G_t     = P * I * P.' * gmm2;
    mpar.I_G_dot_t = zeros(3);
    mpar.M_total   = Mtot * g;

    % =====================================================================
    % 3. THE REST OF THE BASE SEED PARAMS
    % =====================================================================
    bsp.stripGeom = struct('strips', s, 'massParams', mpar, ...
                           'source', sprintf('seedMassModel %s, calibration %s, %d strip(s)/plate', ...
                                             o.name, opts.calibration, nS));
    % planform for the renderers: plates projected on (body z, body x) = (mesh x, mesh y)
    ws = warning('off', 'MATLAB:polyshape:repairedBySimplify');
    shp = polyshape();
    for i = 1:nP
        shp = union(shp, polyshape(pl(i).corners(:,1).' * mm, pl(i).corners(:,2).' * mm));
    end
    warning(ws);
    bsp.seedShape     = shp;
    bsp.seedThickness = sum(thick .* area) / sum(area);
    wingMass          = sum(mUse(~isNut)) * g;
    bsp.seedDensity   = wingMass / sum(area);
    bsp.numStrips     = M;
    bsp.tSamples      = 0;
    bsp.nutPos_t      = P * cen(isNut,:).' * mm;
    bsp.nutMass_t     = sum(mUse(isNut)) * g;

    span = sum(dz);
    meta = struct('name', string(o.name), 'seedNum', o.seedNum, 'calibration', string(opts.calibration), ...
        'measuredMass', o.measuredMass*g, 'predictedMass', sum(mPred)*g, 'calibrationFactor', o.calibration, ...
        'massErrPct', 100*(sum(mPred)/o.measuredMass - 1), 'M_total', Mtot*g, ...
        'nutMass', bsp.nutMass_t, 'wingMass', wingMass, 'nutFrac', bsp.nutMass_t/(Mtot*g), ...
        'span', span, 'area', sum(area), 'AR', span^2/sum(area), ...
        'wingLoading', Mtot*g/sum(area), 'twistRangeDeg', rad2deg(max(twist) - min(twist)), ...
        'maxDihedralDeg', rad2deg(max(abs(dihed))), 'confidence', string(o.confidence));
end

function [a, c] = tri(p1, p2, p3)
% area and centroid of a 3D triangle
    a = 0.5 * norm(cross(p2 - p1, p3 - p1));
    c = (p1 + p2 + p3) / 3;
end
