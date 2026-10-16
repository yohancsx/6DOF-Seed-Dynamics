function result = runSingleMode(name, nutPos, q0, omega0, cfg, baseBsp)
% RUNSINGLEMODE  Build a fixed-nut seed, integrate one drop, and package it.
%
% One named flight-mode case: place the nut at nutPos, release from the given
% orientation/spin, integrate with ode45, and return a result struct (plus a
% one-line metric/classifier summary printed to the console).
%
% INPUTS
%   name    : label for the mode (figure titles / console line).
%   nutPos  : 3x1 nut position in body coords [x;y;z] (m).
%   q0      : 4x1 initial orientation quaternion (body->world, scalar-first).
%   omega0  : 3x1 initial body angular velocity (rad/s).
%   cfg     : config struct (geometry, environment/switches, tspan, tolerances,
%             .aero, .metricOpts). Optional: .seedModeThresholds (default
%             defaultSeedModeThresholds()), .spanLength (default: the built
%             seed's span, sum of the strip widths).
%   baseBsp : baseline base-seed-params (seedShape, density, thickness, numStrips).
%
% OUTPUT
%   result : struct with .name, .t, .x, .seedParams, .metrics (computeSeedModeMetrics),
%            .mode (six-part label, classifySeedMode) and .parts (its six parts,
%            'flip/revolution/attitude/path/regularity').

    % --- Build the seed at this nut position ------------------------------
    bsp = baseBsp;
    bsp.tSamples  = cfg.tSamples;
    bsp.nutPos_t  = repmat(nutPos(:), 1, numel(cfg.tSamples));
    bsp.nutMass_t = cfg.nutMass * ones(size(cfg.tSamples));
    sp = buildSeedParams(bsp, cfg);

    % --- Integrate --------------------------------------------------------
    x0      = [zeros(3,1); q0(:); zeros(3,1); omega0(:)];
    odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);
    rhs     = seedRHS(sp);   % @seed6DOFODE (planar) or @seed6DOFODE3D (shape3d)
    [t, x]  = ode45(@(tt,xx) rhs(tt, xx, sp), cfg.tspan, x0, odeOpts);

    % --- Metrics + six-part classification ---------------------------------
    if isfield(cfg, 'seedModeThresholds') && ~isempty(cfg.seedModeThresholds)
        th = cfg.seedModeThresholds;
    else
        th = defaultSeedModeThresholds();
    end
    mo = cfg.metricOpts;
    if ~isfield(mo, 'spanLength')
        if isfield(cfg, 'spanLength'); mo.spanLength = cfg.spanLength;
        else;                          mo.spanLength = sum(sp.strips.dz); end
    end
    mo.flipHysteresisTurns = th.flipHysteresisTurns;
    result.name       = name;
    result.t          = t;
    result.x          = x;
    result.seedParams = sp;
    result.metrics    = computeSeedModeMetrics(t, x, mo);
    [result.mode, p]  = classifySeedMode(result.metrics, th);
    result.parts      = sprintf('%s/%s/%s/%s/%s', p.flip, p.revolution, p.attitude, ...
                                p.path, p.regularity);

    m = result.metrics;
    fprintf(['%-26s  mode=%-44s  descent=%5.2f  vSpin=%6.2f  spanSpin=%6.2f  ' ...
             'cone=%4.0f  glide=%4.2f  conv=%d\n'], ...
            name, [result.mode ' [' result.parts ']'], m.descentSpeed, m.verticalSpinMag, m.spanwiseSpin, ...
            m.coneAngleDeg, m.glideRatio, m.converged);
end
