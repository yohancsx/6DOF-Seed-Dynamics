function sp = buildSeedParams(bsp, cfg)
% BUILDSEEDPARAMS  Turn a base-seed-params struct into a full, ready-to-run
%   seedParams: geometry + mass via setupSeedShape3D, plus the environment
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
%         .shapeModel (optional) must be 'shape3d' if given (the default). The
%         planar model is archived (archive/planar/, git tag planar-final), so
%         'planar' is an error.
%
% The switches the model honours:
%   enableEdgeDrag, enableAddedMassRate, enableAddedMassMoment, enableAddedMass3D,
%   enableLEV, levApplicationPoint
% An override of one of the archived planar model's switches is WARNED and
% skipped, never silently applied -- an inert knob that looks live is worse than
% a missing one.
%
% LEV TUNABLES: cfg.lev may override any of AR, Kp, Ki, Kv, RoCrit,
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
%   sp  : full seedParams struct accepted by seed6DOFODE3D.

    if isfield(cfg, 'shapeModel') && ~strcmpi(char(cfg.shapeModel), 'shape3d')
        error('buildSeedParams:planarArchived', ...
              ['cfg.shapeModel = ''%s'': the planar model is archived (archive/planar/, ' ...
               'git tag planar-final). Use ''shape3d'' (the default).'], char(cfg.shapeModel));
    end
    sp = setupSeedShape3D(struct('baseSeedParams', bsp));

    sp.rhoFluid = cfg.rhoFluid;
    sp.g        = cfg.g;

    % Physics switches: the BUILDER (setupSeedShape3D) owns the defaults. Only
    % apply EXPLICIT cfg overrides here; never silently clobber a builder default.
    switches = {'enableAddedMass3D', 'enableAddedMassRate', 'enableAddedMassMoment', ...
                'enableEdgeDrag', 'enableLEV', 'levApplicationPoint'};
    for k = 1:numel(switches)
        if isfield(cfg, switches{k})
            sp.(switches{k}) = cfg.(switches{k});
        end
    end

    % The archived planar model's switches: the terms they gated (span force and
    % torque, Tx, Ty) are not in this model. Warn rather than silently ignore, so
    % a config carried over from the planar suites is visible.
    deadSwitches = {'enableSpanTorque', 'enableSpanCOPMigration', ...
                    'enableSpanTorqueAttenuation', 'enableTxDamping', ...
                    'enableNormalSpinDamping', 'enableSpanForce', ...
                    'enableSpanGeomVelocity'};
    for k = 1:numel(deadSwitches)
        if isfield(cfg, deadSwitches{k})
            warning('buildSeedParams:deadSwitch', ...
                ['cfg.%s is ignored: the term it gated belonged to the archived planar ' ...
                 'model; drop it from this config.'], deadSwitches{k});
        end
    end

    % The application point is a string switch: reject a typo here, at build
    % time, rather than deep inside an ode45 call.
    if isfield(sp, 'levApplicationPoint')
        sp.levApplicationPoint = char(sp.levApplicationPoint);
        if ~any(strcmp(sp.levApplicationPoint, {'colocated', 'forward'}))
            error('buildSeedParams:badLevApplicationPoint', ...
                  'levApplicationPoint must be ''colocated'' or ''forward'' (got ''%s'').', ...
                  sp.levApplicationPoint);
        end
    end

    % LEV tunables (see header for the precedence rules).
    if isfield(cfg, 'lev') && ~isempty(cfg.lev)
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

    % The Rossby definition is a string switch, like levApplicationPoint: reject a
    % typo at build time rather than on the first RHS evaluation.
    if isfield(sp, 'lev') && isfield(sp.lev, 'rossbyDefinition')
        sp.lev.rossbyDefinition = char(sp.lev.rossbyDefinition);
        if ~any(strcmp(sp.lev.rossbyDefinition, {'kinematic', 'geometric'}))
            error('buildSeedParams:badRossbyDefinition', ...
                  'lev.rossbyDefinition must be ''kinematic'' or ''geometric'' (got ''%s'').', ...
                  sp.lev.rossbyDefinition);
        end
    end

    % Same for the aero constants those terms carried (now archived with them,
    % in archive/planar/aero/planarAeroConstants.m).
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        deadAero = {'C_span_torque', 'k0_spanTorque', 'C_Tx', 'C_fy', 'C_span'};
        stale    = deadAero(isfield(cfg.aero, deadAero));
        if ~isempty(stale)
            warning('buildSeedParams:deadAero', ...
                ['cfg.aero.%s is ignored: the term it scaled belonged to the archived ' ...
                 'planar model; drop it from this config.'], strjoin(stale, ', cfg.aero.'));
        end
    end

    % Aero coefficient overrides are optional; if absent, computeAeroCoeffs uses
    % its built-in (minimal_imp) defaults.
    if isfield(cfg, 'aero') && ~isempty(cfg.aero)
        sp.aero = cfg.aero;
    end
end
