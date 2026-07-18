function result = runOneSeedCase(cs, baselineSeedParams, baselineBsp, cfg)
% RUNONESEEDCASE  Build, integrate, and package one test-suite case.
%
% Decides whether the case needs a seed rebuild (mass/nut-position/strip-count
% changes) or can reuse the cached baseline seed (initial-condition-only cases),
% applies any post-build multiplier mutation, integrates with ode45, and returns
% a self-contained result struct. Never throws: integration failures are caught
% and recorded in .status so one bad case can't abort the suite.
%
% INPUTS
%   cs                 : one case struct from seedTestCases.
%   baselineSeedParams : the pre-built baseline seedParams (reused by IC-only cases).
%   baselineBsp        : the baseline base-seed-params (merged with overrides on rebuild).
%   cfg                : suite config (tspan, tolerances, environment, ...).
%
% OUTPUT
%   result : struct with fields (fixed order, so results concatenate cleanly):
%     .group .label .sweepVal   - identity/labelling
%     .status                   - 'OK' | 'NONFINITE' | 'FAILED'
%     .errMsg                   - error text if FAILED, else ''
%     .x0                       - 13x1 initial state used
%     .t, .x                    - ode45 time (Nx1) and state (Nx13) history
%     .descended                - true if the seed lost net altitude (world Y)
%     .seedParams               - the exact seedParams integrated (for re-plotting
%                                 / later intermediate recovery)

% --- Initialise all fields up front (fixed order for clean concatenation) ----
result.group      = cs.group;
result.label      = cs.label;
result.sweepVal   = cs.sweepVal;
result.status     = 'OK';
result.errMsg     = '';
result.x0         = [];
result.t          = [];
result.x          = [];
result.descended  = false;
result.seedParams = struct();

% --- Build (or reuse) the seed --------------------------------------------
if isempty(cs.bspOverrides)
    sp = baselineSeedParams;                        % IC-only case: reuse cached seed
else
    bsp = mergeFields(baselineBsp, cs.bspOverrides);% seed changed: rebuild
    sp  = buildSeedParams(bsp, cfg);
end
if ~isempty(cs.postBuildFn)
    sp = cs.postBuildFn(sp);                        % e.g. per-strip lift/drag multipliers
end
result.seedParams = sp;

% --- Initial state x = [r; q; v; omega] -----------------------------------
x0 = [zeros(3,1); cs.q0(:); cs.v0(:); cs.omega0(:)];
result.x0 = x0;

% --- Integrate ------------------------------------------------------------
odeOpts = odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol);
try
    [tOut, xOut] = ode45(@(t, x) seed6DOFODE(t, x, sp), cfg.tspan, x0, odeOpts);
    result.t = tOut;
    result.x = xOut;
    if any(~isfinite(xOut(:)))
        result.status = 'NONFINITE';                % blew up mid-integration
    end
catch ME
    result.status = 'FAILED';
    result.errMsg = ME.message;
end

% --- Descent sanity flag (world Y is state component 2) -------------------
if ~isempty(result.x)
    result.descended = result.x(end, 2) < result.x(1, 2);
end
end


% =========================================================================
% LOCAL: overlay override fields onto a base struct
% =========================================================================
function s = mergeFields(s, ov)
    f = fieldnames(ov);
    for i = 1:numel(f)
        s.(f{i}) = ov.(f{i});
    end
end
