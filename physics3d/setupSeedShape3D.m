function seedParamsFull = setupSeedShape3D(seedParamsIn)
% SETUPSEEDSHAPE3D  Build a NON-PLANAR (3D-shape) seed's geometry + mass.
%
% The 'shape3d' analogue of setupSeedShapeAndMass. It produces the same strip +
% mass-property fields the ODE consumes, PLUS the per-strip 3D contract fields
% (ygc_body and the local chord/normal/span frames) that let seed6DOFODE3D place
% and orient each strip out of the body x-z plane (twist / camber / dihedral).
%
% STATUS: supports flat plates and TWIST (spanwise geometric pitch). It wraps
% setupSeedShapeAndMass for the planar strips + mass (byte-identical to the planar
% builder), then rotates each strip's aerodynamic frame about the span axis by the
% requested twist theta(z) and retags the model 'shape3d'. With no twist the
% frames are identity and the seed is a strict superset of the planar seed (the
% flat-equivalence regression). Out-of-plane CURVATURE / dihedral is not yet
% implemented (see the SHAPE HOOK in the body).
%
% TWIST INPUT: set baseSeedParams.twist to a scalar (uniform pitch, rad), a
% length-M vector (per strip), or a function handle @(z) giving the pitch at
% spanwise position z (m) -- e.g. an anti-symmetric linear twist
% @(z) k*z/(spanLength/2). Absent/empty -> flat.
%
% INPUT
%   seedParamsIn : struct with .baseSeedParams (same as setupSeedShapeAndMass;
%                  a flat polyshape planform + nut mass/position over time, plus
%                  the optional .twist above).
% OUTPUT
%   seedParamsFull : the planar build, plus
%                    .strips.ygc_body  (1xM, = 0; strips stay at y=0 under twist)
%                    .strips.chordDir  (3xM, strip chord  axis in body coords)
%                    .strips.normalDir (3xM, strip normal axis in body coords)
%                    .strips.spanDir   (3xM, strip span   axis in body coords)
%                    .strips.twist     (1xM, the applied pitch per strip, rad)
%                    .model = 'shape3d'

    bsp = seedParamsIn.baseSeedParams;

    % --- Planar base: identical strips + mass to the flat-plate builder -----
    seedParamsFull = setupSeedShapeAndMass(seedParamsIn);

    M   = numel(seedParamsFull.strips.chord);
    zgc = seedParamsFull.strips.zgc_body;        % 1xM strip spanwise centres (body z)

    % --- TWIST: per-strip geometric pitch theta(z) about the LOCAL SPAN axis --
    % Twisting rotates each strip's aerodynamic frame about the span (body z) axis
    % by theta, changing its incidence, while KEEPING the strip centre at y = 0
    % (the chord rotates about the strip's own mid-chord). So the mass/inertia
    % stay the flat-plate values -- twist is treated as an aerodynamic
    % re-orientation of a thin plate, exact in the thin-plate limit -- and ONLY
    % the per-strip frames change. An ANTI-symmetric theta(z) turns some of each
    % strip's lift into a chordwise force at a spanwise arm, producing a net
    % torque about the normal axis: a centred-CoM seed then spins up (the samara
    % autorotation mechanism). theta = 0 gives identity frames (the flat plate).
    %
    % SHAPE HOOK (still TODO): out-of-plane CURVATURE / dihedral -- strips leaving
    % y = 0 -- additionally needs ygc_body != 0, a strip velocity that uses ygc,
    % rotated local inertia, and a full-tensor added mass. Twist needs none of
    % those, which is why it lands here first.
    theta = resolveTwist(bsp, zgc, M);           % 1xM pitch per strip (rad); 0 -> flat
    c = cos(theta);   s = sin(theta);            % 1xM

    seedParamsFull.strips.ygc_body  = zeros(1, M);               % strips remain at body y = 0
    seedParamsFull.strips.chordDir  = [ c;  s; zeros(1, M)];     % Rz(theta) * [1;0;0]
    seedParamsFull.strips.normalDir = [-s;  c; zeros(1, M)];     % Rz(theta) * [0;1;0]
    seedParamsFull.strips.spanDir   = [zeros(2, M); ones(1, M)]; % Rz(theta) * [0;0;1] = span
    seedParamsFull.strips.twist     = theta;                     % record the twist profile (rad)

    % --- Retag the model (setupSeedShapeAndMass stamped it 'planar') --------
    seedParamsFull.model = 'shape3d';
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
