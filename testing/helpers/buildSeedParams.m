function sp = buildSeedParams(bsp, cfg)
% BUILDSEEDPARAMS  Turn a base-seed-params struct into a full, ready-to-run
%   seedParams: geometry + mass via setupSeedShapeAndMass, plus the environment
%   fields seed6DOFODE needs (rhoFluid, g, enableSpanForce, optional aero).
%
% Shared by the test-suite driver (for the baseline seed) and runOneSeedCase
% (for cases that modify the seed and must rebuild).
%
% INPUTS
%   bsp : base seed params (the .baseSeedParams sub-struct for setupSeedShapeAndMass).
%   cfg : suite config struct; uses .rhoFluid, .g, and the optional
%         .enableSpanForce / .enableSpanGeomVelocity (both DEFAULT TRUE if
%         absent) and .aero (coefficient overrides for computeAeroCoeffs).
%         .shapeModel (optional): 'planar' (default) or 'shape3d'. 'shape3d'
%         routes to setupSeedShape3D (needs physics3d/ on the path) and produces
%         the non-planar seed contract; 'planar' uses setupSeedShapeAndMass.
%
% OUTPUT
%   sp  : full seedParams struct accepted by seed6DOFODE (planar) or
%         seed6DOFODE3D (shape3d).

    if isfield(cfg, 'shapeModel') && strcmpi(char(cfg.shapeModel), 'shape3d')
        sp = setupSeedShape3D(struct('baseSeedParams', bsp));
    else
        sp = setupSeedShapeAndMass(struct('baseSeedParams', bsp));
    end

    sp.rhoFluid = cfg.rhoFluid;
    sp.g        = cfg.g;

    % Physics switches: the BUILDER owns the defaults -- setupSeedShapeAndMass
    % stamps the "FULL-minus-geomVelocity" config, and setupSeedShape3D adjusts it
    % for the shape (e.g. it turns the planar span-force hack OFF for a curved
    % seed, where the tilted strips already supply that force from geometry). So
    % only apply EXPLICIT cfg overrides here; never silently clobber a
    % builder-chosen default. (For a planar seed with no overrides this leaves
    % exactly the same values this function used to hardcode.)
    switches = {'enableSpanForce', 'enableSpanTorque', 'enableSpanGeomVelocity', ...
                'enableSpanCOPMigration', 'enableSpanTorqueAttenuation', ...
                'enableTxDamping', 'enableAddedMass3D'};
    for k = 1:numel(switches)
        if isfield(cfg, switches{k})
            sp.(switches{k}) = cfg.(switches{k});
        end
    end

    % Aero coefficient overrides are optional; if absent, computeAeroCoeffs uses
    % its built-in (minimal_imp) defaults.
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        sp.aero = cfg.aero;
    end
end
