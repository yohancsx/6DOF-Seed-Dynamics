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
%   cfg : suite config struct; uses .rhoFluid, .g, .enableSpanForce, and the
%         optional .aero (empirical coefficient overrides for computeAeroCoeffs).
%
% OUTPUT
%   sp  : full seedParams struct accepted by seed6DOFODE.

    sp = setupSeedShapeAndMass(struct('baseSeedParams', bsp));

    sp.rhoFluid        = cfg.rhoFluid;
    sp.g               = cfg.g;
    sp.enableSpanForce = cfg.enableSpanForce;

    % Aero coefficient overrides are optional; if absent, computeAeroCoeffs uses
    % its built-in (minimal_imp) defaults.
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        sp.aero = cfg.aero;
    end
end
