function info = validateSeedParams(seedParams)
% VALIDATESEEDPARAMS  Check a seedParams struct against the seed-model contract.
%
% The contract is the interface between a seed-GEOMETRY model (which builds the
% strips + mass properties for a particular shape) and the shape-AGNOSTIC code
% that consumes them: the ODE right-hand side, getMassProperties / getAddedMass,
% computeSeedLocalVel, and the rigid-body core (rigidBody6DOF). Both the planar
% flat-plate builder (setupSeedShapeAndMass -> model 'planar') and the future
% full-3D builder (-> model 'shape3d') must satisfy it, so the shared code can
% rely on a fixed set of fields regardless of shape.
%
% This is a DIAGNOSTIC / test helper -- it is NOT called inside the ODE RHS (that
% runs thousands of times per integration). Call it once after building a seed,
% or in tests, to catch a malformed struct early.
%
% =========================================================================
% THE CONTRACT
% -------------------------------------------------------------------------
% seedParams.model                : 'planar' | 'shape3d'   (the model tag)
% seedParams.rhoFluid             : scalar, fluid density (kg/m^3)
% seedParams.g                    : scalar, gravity magnitude (m/s^2)
% seedParams.strips.chord         : 1xM strip chord lengths (m)
% seedParams.strips.dz            : 1xM strip WIDTHS (m). Planar: spanwise (body z)
%                                   extent; 3D: arc-length of the strip along the
%                                   (curved) span -- i.e. the true strip width, so
%                                   area = chord.*dz is correct for both.
% seedParams.strips.xgc_body      : 1xM strip geometric-centre, body x (chordwise)
% seedParams.strips.zgc_body      : 1xM strip geometric-centre, body z (spanwise)
% seedParams.strips.liftMult      : 1xM per-strip translational-lift multipliers
% seedParams.strips.dragMult      : 1xM per-strip drag multipliers
% seedParams.massParams.tSamples  : Nx1 time samples for the mass properties (s)
% seedParams.massParams.com_t     : 3xN CoM in body-datum coords vs time (m)
% seedParams.massParams.I_G_t     : 3x3xN inertia about the CoM vs time (kg m^2)
% seedParams.massParams.I_G_dot_t : 3x3xN d/dt of I_G vs time (kg m^2/s)
% seedParams.massParams.M_total   : scalar total mass (kg)
%
% MODEL-SPECIFIC:
%   'planar' : nothing extra. Strips lie in the body x-z plane (y = 0) and the
%              aerodynamic normal is body y; the RHS hardcodes those.
%   'shape3d': additionally, each strip carries its 3D position and its local
%              orthonormal frame in BODY coordinates, so the RHS can project the
%              strip velocity into the strip's own (chord, normal) plane and
%              rotate the resulting force back to body:
%                seedParams.strips.ygc_body  : 1xM strip geometric-centre, body y
%                seedParams.strips.chordDir  : 3xM strip chord  axis, body coords
%                seedParams.strips.normalDir : 3xM strip normal axis, body coords
%                seedParams.strips.spanDir   : 3xM strip span   axis, body coords
%              (Planar is the special case chordDir=[1;0;0], normalDir=[0;1;0],
%              spanDir=[0;0;1], ygc_body=0 for every strip.)
%
% OPTIONAL (defaults applied by the RHS if absent): seedParams.aero, and the
% physics switches enableSpanForce / enableSpanTorque / enableSpanGeomVelocity /
% enableSpanCOPMigration / enableSpanTorqueAttenuation / enableTxDamping.
%
% INPUT
%   seedParams : struct to check.
% OUTPUT
%   info : struct with .model, .numStrips, .ok (true if no error was thrown).
% Throws an error (id 'validateSeedParams:*') on the first violation.

    mustBeField(seedParams, 'model', 'seedParams');
    model = seedParams.model;
    if ~ischar(model) && ~isstring(model)
        error('validateSeedParams:badModel', '.model must be a char/string, got %s.', class(model));
    end
    model = char(model);

    % --- Base scalars -----------------------------------------------------
    mustBeField(seedParams, 'rhoFluid', 'seedParams');
    mustBeField(seedParams, 'g',        'seedParams');

    % --- Strip geometry ---------------------------------------------------
    mustBeField(seedParams, 'strips', 'seedParams');
    s = seedParams.strips;
    stripFields = {'chord','dz','xgc_body','zgc_body','liftMult','dragMult'};
    for f = stripFields
        mustBeField(s, f{1}, 'seedParams.strips');
    end
    M = numel(s.chord);
    for f = stripFields
        if numel(s.(f{1})) ~= M
            error('validateSeedParams:stripLen', ...
                  'strips.%s has length %d, expected %d (= numel(strips.chord)).', ...
                  f{1}, numel(s.(f{1})), M);
        end
    end

    % --- Mass properties --------------------------------------------------
    mustBeField(seedParams, 'massParams', 'seedParams');
    mp = seedParams.massParams;
    for f = {'tSamples','com_t','I_G_t','I_G_dot_t','M_total'}
        mustBeField(mp, f{1}, 'seedParams.massParams');
    end
    N = numel(mp.tSamples);
    assertSize(mp.com_t,     [3 N],     'massParams.com_t');
    assertSize(mp.I_G_t,     [3 3 N],   'massParams.I_G_t');
    assertSize(mp.I_G_dot_t, [3 3 N],   'massParams.I_G_dot_t');

    % --- Model-specific ---------------------------------------------------
    switch model
        case 'planar'
            % nothing extra
        case 'shape3d'
            extra = {'ygc_body','chordDir','normalDir','spanDir'};
            for f = extra
                mustBeField(s, f{1}, 'seedParams.strips (shape3d)');
            end
            if numel(s.ygc_body) ~= M
                error('validateSeedParams:stripLen', 'strips.ygc_body has length %d, expected %d.', numel(s.ygc_body), M);
            end
            for f = {'chordDir','normalDir','spanDir'}
                assertSize(s.(f{1}), [3 M], sprintf('strips.%s', f{1}));
            end
        otherwise
            error('validateSeedParams:unknownModel', ...
                  'Unknown seedParams.model "%s" (expected ''planar'' or ''shape3d'').', model);
    end

    info = struct('model', model, 'numStrips', M, 'ok', true);
end


% =========================================================================
% LOCAL helpers
% =========================================================================
function mustBeField(s, name, ctx)
    if ~isstruct(s) || ~isfield(s, name)
        error('validateSeedParams:missingField', 'Missing required field %s.%s', ctx, name);
    end
end

function assertSize(v, sz, name)
    actual = size(v, 1:numel(sz));
    if ~isequal(actual, sz)
        error('validateSeedParams:badSize', '%s has size [%s], expected [%s].', ...
              name, num2str(actual), num2str(sz));
    end
end
