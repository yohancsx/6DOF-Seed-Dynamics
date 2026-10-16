function sp = buildPlanarSeedParams(bsp, cfg)
% BUILDPLANARSEEDPARAMS  Build an ARCHIVED planar seed (for comparison tests only).
%
% The planar branch of testing/helpers/buildSeedParams as it was before the planar
% model was archived (2026-10-15; the original is at git tag planar-final): the
% planar builder setupSeedShapeAndMass, the environment fields, and explicit
% overrides of the planar model's own switches. Its RHS is seed6DOFODE, in this
% same folder.
%
% Only regression tests call this, to keep checking that the active 3D model
% reproduces the planar model's CORE exactly (testShape3DFlatEquivalence and the
% planar halves of testCurvedMassProperties / testNewPhysicsTerms). They put
% archive/planar/physics and archive/planar/aero on the path themselves.
%
% INPUTS  bsp : base seed params (the planar builder's .baseSeedParams).
%         cfg : .rhoFluid, .g, optional planar switch overrides, optional .aero.
% OUTPUT  sp  : planar seedParams (model 'planar'), for seed6DOFODE.

    sp = setupSeedShapeAndMass(struct('baseSeedParams', bsp));
    sp.rhoFluid = cfg.rhoFluid;
    sp.g        = cfg.g;
    switches = {'enableSpanForce', 'enableSpanTorque', 'enableSpanGeomVelocity', ...
                'enableSpanCOPMigration', 'enableSpanTorqueAttenuation', ...
                'enableTxDamping', 'enableNormalSpinDamping', 'enableAddedMass3D', ...
                'enableAddedMassRate', 'enableAddedMassMoment'};
    for k = 1:numel(switches)
        if isfield(cfg, switches{k}); sp.(switches{k}) = cfg.(switches{k}); end
    end
    if isfield(cfg, 'aero') && ~isempty(cfg.aero); sp.aero = cfg.aero; end
end
