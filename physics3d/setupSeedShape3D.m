function seedParamsFull = setupSeedShape3D(seedParamsIn)
% SETUPSEEDSHAPE3D  Build a NON-PLANAR (3D-shape) seed's geometry + mass.
%
% The 'shape3d' analogue of setupSeedShapeAndMass. It produces the same strip +
% mass-property fields the ODE consumes, PLUS the per-strip 3D contract fields
% (ygc_body and the local chord/normal/span frames) that let seed6DOFODE3D place
% and orient each strip out of the body x-z plane (twist / camber / dihedral).
%
% STATUS: supports flat plates, TWIST (spanwise geometric pitch about the span
% axis), and spanwise CURVATURE (dihedral bending about the chord axis, strips
% leaving y=0). It wraps setupSeedShapeAndMass for the planar strips, then sets
% each strip's out-of-plane position + local frame from the requested twist and
% curvature and retags 'shape3d'. With neither, the frames are identity and the
% seed is a strict superset of the planar seed (the flat-equivalence regression).
%
% IMPORTANT (mass/inertia not yet curved -- step 2): the mass properties come
% straight from setupSeedShapeAndMass (a FLAT plate), so for a CURVED seed the
% CoM and inertia are still the flat-plate values and its DYNAMICS are NOT
% physically correct yet. Its geometry, per-strip velocity, and rendering ARE
% correct. Twist keeps strips at y=0, so a twisted seed's mass is exact.
%
% TWIST INPUT: baseSeedParams.twist = scalar (uniform pitch, rad) / length-M
% vector / handle @(z) (pitch vs spanwise position) -- e.g. anti-symmetric
% @(z) k*z. Absent -> no twist.
% CURVATURE INPUT: baseSeedParams.curvature = the dihedral tangent angle phi vs
% ARC LENGTH s, same forms -- e.g. @(s) k*s for a circular-arc bowl, or a scalar
% for a uniform tilt. Absent -> flat (no bending). Twist and curvature compose.
%
% INPUT
%   seedParamsIn : struct with .baseSeedParams (a flat polyshape planform + nut
%                  mass/position over time, plus optional .twist / .curvature).
% OUTPUT
%   seedParamsFull : the planar build, plus
%                    .strips.ygc_body  (1xM, out-of-plane strip position, body y)
%                    .strips.chordDir/normalDir/spanDir (3xM strip frame, body)
%                    .strips.twist     (1xM applied pitch per strip, rad)
%                    .strips.dihedral  (1xM applied tangent angle per strip, rad)
%                    .model = 'shape3d'
%                  (.strips.zgc_body/z_body are remapped to the physical spanwise
%                  position on the curve when curvature is set.)

    bsp = seedParamsIn.baseSeedParams;

    % --- Planar base: identical strips + mass to the flat-plate builder -----
    seedParamsFull = setupSeedShapeAndMass(seedParamsIn);

    M    = numel(seedParamsFull.strips.chord);
    sArc = seedParamsFull.strips.zgc_body;       % flat span coordinate = ARC LENGTH
                                                  % along the (bent) profile (1xM)

    % --- TWIST: per-strip geometric pitch theta(s) about the LOCAL SPAN axis ---
    % Twist rotates each strip's frame about the span axis by theta, changing its
    % incidence while KEEPING its centre in place. An anti-symmetric theta(s)
    % turns some of each strip's lift into a chordwise force at a spanwise arm ->
    % a centred-CoM seed spins up (the samara autorotation mechanism).
    theta = resolveTwist(bsp, sArc, M);          % 1xM pitch per strip (rad); 0 -> flat

    % --- CURVATURE: spanwise dihedral tangent angle phi(s) about the CHORD axis -
    % The flat span coordinate sArc is treated as ARC LENGTH; a curvature bends the
    % profile in the body (y,z) plane. phi(s) is the local tangent (dihedral) angle;
    % integrating dz=cos(phi)ds, dy=sin(phi)ds gives each strip's PHYSICAL (z,y)
    % centre on the curve, anchored so s=0 -> (0,0). Strip widths (dz) are the
    % arc-length segments, so they are UNCHANGED by bending. Absent -> flat
    % (phi=0, y=0, positions unchanged), so twist-only and planar stay bit-identical.
    %
    % NOTE (step 2 pending): the mass properties (com_t, I_G_t) are still the
    % FLAT-plate values from setupSeedShapeAndMass. Curvature-correct mass/inertia
    % (the wing-centroid y shift + rotated local inertia) is NOT yet implemented,
    % so a CURVED seed's DYNAMICS are not physically correct until then; its
    % geometry, per-strip velocity, and rendering ARE correct.
    if isfield(bsp, 'curvature') && ~isempty(bsp.curvature)
        [phi, zCurve, yCurve] = curveFromProfile(bsp, sArc, M);
        seedParamsFull.strips.zgc_body = zCurve;    % physical spanwise position (body z)
        seedParamsFull.strips.z_body   = zCurve;
        ygc = yCurve;
    else
        phi = zeros(1, M);
        ygc = zeros(1, M);
    end

    % --- Per-strip frame = dihedral(phi, about chord) o twist(theta, about span).
    % Columns of R_strip = Rd(phi)*Rz(theta), in body coords. Reduces to the twist
    % frame when phi=0 and to identity when theta=phi=0 (so flat/twist unchanged).
    ct = cos(theta);   st = sin(theta);   cp = cos(phi);   sp = sin(phi);   % 1xM
    seedParamsFull.strips.chordDir  = [ct;         cp.*st;   -sp.*st];
    seedParamsFull.strips.normalDir = [-st;        cp.*ct;   -sp.*ct];
    seedParamsFull.strips.spanDir   = [zeros(1,M); sp;        cp];
    seedParamsFull.strips.ygc_body  = ygc;
    seedParamsFull.strips.twist     = theta;     % record profiles (rad)
    seedParamsFull.strips.dihedral  = phi;

    % --- Retag the model (setupSeedShapeAndMass stamped it 'planar') --------
    seedParamsFull.model = 'shape3d';
end


% =========================================================================
% LOCAL: curvature profile -> per-strip tangent angle + physical (z,y) centres
% =========================================================================
function [phi, zNew, yNew] = curveFromProfile(bsp, sCenters, M)
% bsp.curvature (same forms as bsp.twist) gives the dihedral tangent angle phi as
% a function of ARC LENGTH s: a scalar (uniform tilt -> a shallow V), a length-M
% vector (per strip), or a handle @(s) (e.g. @(s) k*s for a circular-arc bowl).
% Returns phi at the strip centres and the physical (z,y) centre of each strip on
% the integrated curve (arc-length s=0 anchored to the origin).
    cv = bsp.curvature;
    ev = @(sq) evalProfile(cv, sq, sCenters, M);   % phi at query arc-lengths
    phi = ev(sCenters);

    sGrid = linspace(min(sCenters), max(sCenters), 1001);
    phiG  = ev(sGrid);
    zG = cumtrapz(sGrid, cos(phiG));               % dz = cos(phi) ds
    yG = cumtrapz(sGrid, sin(phiG));               % dy = sin(phi) ds
    zG = zG - interp1(sGrid, zG, 0, 'linear', 'extrap');   % anchor s=0 -> (0,0)
    yG = yG - interp1(sGrid, yG, 0, 'linear', 'extrap');
    zNew = interp1(sGrid, zG, sCenters, 'linear', 'extrap');
    yNew = interp1(sGrid, yG, sCenters, 'linear', 'extrap');
end

function v = evalProfile(cv, sq, sCenters, M)
% Evaluate a scalar / length-M vector / handle profile at arc-lengths sq.
    if isa(cv, 'function_handle')
        v = arrayfun(@(s) cv(s), sq(:).');
    elseif isscalar(cv)
        v = double(cv) * ones(1, numel(sq));
    elseif numel(cv) == M
        v = interp1(sCenters, double(cv(:)).', sq, 'linear', 'extrap');
    else
        error('setupSeedShape3D:badCurvature', ...
              'bsp.curvature must be a scalar, a length-%d vector, or a handle @(s).', M);
    end
end


% =========================================================================
% LOCAL: resolve the twist input into a 1xM per-strip pitch (rad)
% =========================================================================
function theta = resolveTwist(bsp, zgc, M)
% bsp.twist may be: absent/empty (flat, theta = 0); a scalar (uniform pitch on
% every strip); a length-M vector (per-strip pitch, strip order); or a function
% handle @(z) giving the pitch at spanwise position z (m), evaluated at each
% strip's spanwise centre -- e.g. a linear anti-symmetric twist
% @(z) k*z/(spanLength/2).
    if ~isfield(bsp, 'twist') || isempty(bsp.twist)
        theta = zeros(1, M);
        return
    end
    tw = bsp.twist;
    if isa(tw, 'function_handle')
        theta = arrayfun(@(z) tw(z), zgc(:).');          % 1xM
    elseif isscalar(tw)
        theta = double(tw) * ones(1, M);
    elseif numel(tw) == M
        theta = double(tw(:)).';
    else
        error('setupSeedShape3D:badTwist', ...
              'bsp.twist must be a scalar, a length-%d vector, or a function handle @(z).', M);
    end
end
