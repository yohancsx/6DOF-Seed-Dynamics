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
%             enableNormalSpinDamping, enableAddedMass3D, enableAddedMassRate,
%             enableAddedMassMoment
%   shape3d : enableEdgeDrag, enableAddedMassRate, enableAddedMassMoment,
%             enableAddedMass3D, enableLEV, levApplicationPoint
% An override for a switch the chosen model does not honour is WARNED and
% skipped, never silently applied -- an inert knob that looks live is worse than
% a missing one.
%
% LEV TUNABLES (shape3d only): cfg.lev may override any of AR, Kp, Ki, Kv, RoCrit,
% p, lambdaV, rossbyDefinition on the builder-computed seedParams.lev. Precedence,
% so the constants stay consistent with Rezgui et al. (2020) eq. (4),
% Kv = Kp - Kp^2*Ki:
%   1. cfg.lev.AR re-derives Kp, Ki and Kv from the new aspect ratio;
%   2. explicit cfg.lev.Kp / .Ki / .Kv values then win;
%   3. if Kp or Ki was set but Kv was not, Kv is recomputed from eq. (4).
% A resulting Kv <= 0 (possible only with inconsistent Kp and Ki) is warned.
% cfg.lev.rossbyDefinition ('kinematic' default | 'geometric') selects which Rossby
% number the gate reads; like levApplicationPoint it is validated here, at build
% time, rather than deep inside an ode45 call.
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
                'enableAddedMassRate', 'enableAddedMassMoment', 'enableEdgeDrag', ...
                'enableLEV', 'levApplicationPoint'};

    % Switches the 3D model no longer honours: the Sep-2026 audit removed the
    % span torque, Tx and Ty; phase 4 retired the span force (measured inert) and
    % its velocity-sampling switch. Warn rather than silently apply an inert
    % override, so a config carried over from the planar suites is visible.
    deadInShape3D = {'enableSpanTorque', 'enableSpanCOPMigration', ...
                     'enableSpanTorqueAttenuation', 'enableTxDamping', ...
                     'enableNormalSpinDamping', 'enableSpanForce', ...
                     'enableSpanGeomVelocity'};
    % ...and the ones the planar model never had (these terms live in physics3d/).
    deadInPlanar  = {'enableEdgeDrag', 'enableLEV', 'levApplicationPoint'};
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

    % The application point is a string switch: reject a typo here, at build
    % time, rather than deep inside an ode45 call.
    if isShape3D && isfield(sp, 'levApplicationPoint')
        sp.levApplicationPoint = char(sp.levApplicationPoint);
        if ~any(strcmp(sp.levApplicationPoint, {'colocated', 'forward'}))
            error('buildSeedParams:badLevApplicationPoint', ...
                  'levApplicationPoint must be ''colocated'' or ''forward'' (got ''%s'').', ...
                  sp.levApplicationPoint);
        end
    end

    % LEV tunables (see header for the precedence rules).
    if isfield(cfg, 'lev') && ~isempty(cfg.lev)
        if ~isShape3D
            warning('buildSeedParams:deadSwitch', ...
                'cfg.lev is ignored by the planar model (LEV exists only in the shape3d model).');
        else
            ov    = cfg.lev;
            known = {'AR', 'Kp', 'Ki', 'Kv', 'RoCrit', 'p', 'lambdaV', 'rossbyDefinition'};
            given = fieldnames(ov);
            bad   = given(~ismember(given, known));
            if ~isempty(bad)
                warning('buildSeedParams:unknownLevField', ...
                    'cfg.lev.%s is not an LEV parameter (known: %s); ignored.', ...
                    strjoin(bad, ', cfg.lev.'), strjoin(known, ', '));
            end
            if isfield(ov, 'AR')                           % 1. re-derive from AR
                base = levPlanformConstants(ov.AR);
                for f = {'AR', 'a0', 'Kp', 'Ki', 'Kv'}
                    sp.lev.(f{1}) = base.(f{1});
                end
            end
            for f = intersect(given, setdiff(known, {'AR'})).'   % 2. explicit values win
                sp.lev.(f{1}) = ov.(f{1});
            end
            if (isfield(ov, 'Kp') || isfield(ov, 'Ki')) && ~isfield(ov, 'Kv')
                sp.lev.Kv = sp.lev.Kp - sp.lev.Kp^2 * sp.lev.Ki;  % 3. keep eq. (4)
            end
            if sp.lev.Kv <= 0
                warning('buildSeedParams:nonPositiveKv', ...
                    ['LEV Kv = %.3g <= 0 (Kp = %.3g, Ki = %.3g): the vortex lift would ' ...
                     'oppose the flow. Kp and Ki must come from the same wing model.'], ...
                    sp.lev.Kv, sp.lev.Kp, sp.lev.Ki);
            end
        end
    end

    % The Rossby definition is a string switch, like levApplicationPoint: reject a
    % typo at build time rather than on the first RHS evaluation.
    if isShape3D && isfield(sp, 'lev') && isfield(sp.lev, 'rossbyDefinition')
        sp.lev.rossbyDefinition = char(sp.lev.rossbyDefinition);
        if ~any(strcmp(sp.lev.rossbyDefinition, {'kinematic', 'geometric'}))
            error('buildSeedParams:badRossbyDefinition', ...
                  'lev.rossbyDefinition must be ''kinematic'' or ''geometric'' (got ''%s'').', ...
                  sp.lev.rossbyDefinition);
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
