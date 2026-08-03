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
%
% OUTPUT
%   sp  : full seedParams struct accepted by seed6DOFODE.

    sp = setupSeedShapeAndMass(struct('baseSeedParams', bsp));

    sp.rhoFluid = cfg.rhoFluid;
    sp.g        = cfg.g;

    % Physics switches: fall back to the "FULL-minus-geomVelocity" default
    % (see setupSeedShapeAndMass) when cfg does not specify them.
    if isfield(cfg, 'enableSpanForce')
        sp.enableSpanForce = cfg.enableSpanForce;
    else
        sp.enableSpanForce = true;
    end
    if isfield(cfg, 'enableSpanTorque')
        sp.enableSpanTorque = cfg.enableSpanTorque;
    else
        sp.enableSpanTorque = true;
    end
    if isfield(cfg, 'enableSpanGeomVelocity')
        sp.enableSpanGeomVelocity = cfg.enableSpanGeomVelocity;
    else
        sp.enableSpanGeomVelocity = false;   % no measurable effect on any mode
    end
    if isfield(cfg, 'enableSpanCOPMigration')
        sp.enableSpanCOPMigration = cfg.enableSpanCOPMigration;
    else
        sp.enableSpanCOPMigration = true;
    end
    if isfield(cfg, 'enableSpanTorqueAttenuation')
        sp.enableSpanTorqueAttenuation = cfg.enableSpanTorqueAttenuation;
    else
        sp.enableSpanTorqueAttenuation = false;
    end
    if isfield(cfg, 'enableTxDamping')
        sp.enableTxDamping = cfg.enableTxDamping;
    else
        sp.enableTxDamping = true;
    end

    % Aero coefficient overrides are optional; if absent, computeAeroCoeffs uses
    % its built-in (minimal_imp) defaults.
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        sp.aero = cfg.aero;
    end
end
