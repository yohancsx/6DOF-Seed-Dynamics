function sp = buildSeedParams(bsp, cfg)
% BUILDSEEDPARAMS  Turn a base-seed-params struct into a full, ready-to-run
%   seedParams: geometry + mass via the model's builder, plus the environment
%   fields the RHS needs (rhoFluid, g, optional aero) and any explicit switch
%   overrides.
%
% Shared by the test-suite driver (for the baseline seed) and runOneSeedCase
% (for cases that modify the seed and must rebuild).
%
% INPUTS
%   bsp : base seed params (the .baseSeedParams sub-struct for setupSeedShapeAndMass).
%   cfg : suite config struct; uses .rhoFluid, .g, optional explicit enable*
%         switch overrides, and .aero (coefficient overrides for computeAeroCoeffs).
%         .shapeModel (optional): 'planar' (default) or 'shape3d'. 'shape3d'
%         routes to setupSeedShape3D (needs physics3d/ on the path) and produces
%         the non-planar seed contract; 'planar' uses setupSeedShapeAndMass.
%
% The two models honour DIFFERENT switch sets:
%   planar  : enableSpanForce, enableSpanGeomVelocity, enableSpanTorque,
%             enableSpanCOPMigration, enableSpanTorqueAttenuation, enableTxDamping,
%             enableNormalSpinDamping, enableAddedMass3D, enableAddedMassRate
%   shape3d : enableEdgeDrag, enableAddedMassRate, enableAddedMass3D
% An override for a switch the chosen model does not honour is WARNED and
% skipped, never silently applied -- an inert knob that looks live is worse than
% a missing one.
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
    % stamps the planar config and setupSeedShape3D stamps the shape3d one. So
    % only apply EXPLICIT cfg overrides here; never silently clobber a
    % builder-chosen default. (For a planar seed with no overrides this leaves
    % exactly the same values this function used to hardcode.)
    switches = {'enableSpanForce', 'enableSpanTorque', 'enableSpanGeomVelocity', ...
                'enableSpanCOPMigration', 'enableSpanTorqueAttenuation', ...
                'enableTxDamping', 'enableNormalSpinDamping', 'enableAddedMass3D', ...
                'enableAddedMassRate', 'enableEdgeDrag'};

    % Switches the 3D model no longer honours: the Sep-2026 audit removed the
    % span torque, Tx and Ty; phase 4 retired the span force (measured inert) and
    % its velocity-sampling switch. Warn rather than silently apply an inert
    % override, so a config carried over from the planar suites is visible.
    deadInShape3D = {'enableSpanTorque', 'enableSpanCOPMigration', ...
                     'enableSpanTorqueAttenuation', 'enableTxDamping', ...
                     'enableNormalSpinDamping', 'enableSpanForce', ...
                     'enableSpanGeomVelocity'};
    % ...and the one the planar model never had (edge drag lives in physics3d/).
    deadInPlanar  = {'enableEdgeDrag'};
    isShape3D = strcmpi(sp.model, 'shape3d');

    for k = 1:numel(switches)
        if isfield(cfg, switches{k})
            if isShape3D && any(strcmp(switches{k}, deadInShape3D))
                warning('buildSeedParams:deadSwitch', ...
                    ['cfg.%s is ignored by the shape3d model (the term it gated was ' ...
                     'removed); drop it from this config.'], switches{k});
                continue
            end
            if ~isShape3D && any(strcmp(switches{k}, deadInPlanar))
                warning('buildSeedParams:deadSwitch', ...
                    ['cfg.%s is ignored by the planar model (the term exists only in ' ...
                     'the shape3d model); drop it or set cfg.shapeModel = ''shape3d''.'], ...
                    switches{k});
                continue
            end
            sp.(switches{k}) = cfg.(switches{k});
        end
    end

    % Same for the aero constants those terms carried. They still exist in the
    % shared computeAeroCoeffs because the frozen planar model reads them.
    if isShape3D && isfield(cfg, 'aero') && ~isempty(cfg.aero)
        deadAero = {'C_span_torque', 'k0_spanTorque', 'C_Tx', 'C_fy', 'C_span'};
        stale    = deadAero(isfield(cfg.aero, deadAero));
        if ~isempty(stale)
            warning('buildSeedParams:deadAero', ...
                ['cfg.aero.%s is ignored by the shape3d model (the term it scaled was ' ...
                 'removed); drop it from this config.'], strjoin(stale, ', cfg.aero.'));
        end
    end

    % Aero coefficient overrides are optional; if absent, computeAeroCoeffs uses
    % its built-in (minimal_imp) defaults.
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        sp.aero = cfg.aero;
    end
end
