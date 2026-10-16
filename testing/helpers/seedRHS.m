function rhs = seedRHS(seedParams)
% SEEDRHS  The ODE right-hand-side function for a seed.
%
% Returns @seed6DOFODE3D for a 'shape3d' seed (the only active model), so the
% integrating helpers and drivers have one place that picks the RHS. Build the
% seed with buildSeedParams (which stamps .model), then integrate with:
%     rhs = seedRHS(sp);
%     [t, x] = ode45(@(t, x) rhs(t, x, sp), tspan, x0, opts);
%
% A 'planar' seed is refused: the planar model was archived on 2026-10-15
% (archive/planar/, git tag planar-final). A struct without .model is refused
% too -- it can only have come from the archived planar builder.
%
% INPUT   seedParams : struct with .model = 'shape3d'.
% OUTPUT  rhs        : function handle f(t, x, seedParams) -> xdot.

    if isfield(seedParams, 'model')
        model = char(seedParams.model);
    else
        model = '(none)';
    end

    switch model
        case 'shape3d'
            rhs = @seed6DOFODE3D;
        case {'planar', '(none)'}
            error('seedRHS:planarArchived', ...
                  ['A planar seed (model "%s"): the planar model is archived ' ...
                   '(archive/planar/, git tag planar-final). Build with buildSeedParams.'], model);
        otherwise
            error('seedRHS:unknownModel', ...
                  'Unknown seedParams.model "%s" (expected ''shape3d'').', model);
    end
end
