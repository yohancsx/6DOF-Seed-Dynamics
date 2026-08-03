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
%             .aero, .metricOpts, .modeThresholds).
%   baseBsp : baseline base-seed-params (seedShape, density, thickness, numStrips).
%
% OUTPUT
%   result : struct with .name, .t, .x, .seedParams, .metrics, .mode.

    % --- Build the seed at this nut position ------------------------------
    bsp = baseBsp;
    bsp.tSamples  = cfg.tSamples;
    bsp.nutPos_t  = repmat(nutPos(:), 1, numel(cfg.tSamples));
    bsp.nutMass_t = cfg.nutMass * ones(size(cfg.tSamples));
    sp = buildSeedParams(bsp, cfg);

    % --- Integrate --------------------------------------------------------
    x0      = [zeros(3,1); q0(:); zeros(3,1); omega0(:)];
    odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);
    [t, x]  = ode45(@(tt,xx) seed6DOFODE(tt, xx, sp), cfg.tspan, x0, odeOpts);

    % --- Metrics + (reference) classification -----------------------------
    result.name       = name;
    result.t          = t;
    result.x          = x;
    result.seedParams = sp;
    result.metrics    = computeTrajectoryMetrics(t, x, cfg.metricOpts);
    [result.mode, ~]  = classifyFlightMode(result.metrics, cfg.modeThresholds);

    m = result.metrics;
    fprintf(['%-26s  guess=%-14s  descent=%5.2f  vSpin=%6.2f  spanSpin=%6.2f  ' ...
             'cone=%4.0f  glide=%4.2f  conv=%d\n'], ...
            name, result.mode, m.descentSpeed, m.verticalSpinMag, m.spanwiseSpin, ...
            m.coneAngleDeg, m.glideRatio, m.converged);
end
