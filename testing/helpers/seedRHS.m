function rhs = seedRHS(seedParams)
% SEEDRHS  The ODE right-hand-side function for a seed, chosen by its model tag.
%
% Returns @seed6DOFODE for a 'planar' seed and @seed6DOFODE3D for a 'shape3d'
% seed, so the integrating helpers and drivers stay model-agnostic. Build the
% seed with buildSeedParams (which stamps .model), then integrate with:
%     rhs = seedRHS(sp);
%     [t, x] = ode45(@(t, x) rhs(t, x, sp), tspan, x0, opts);
%
% The returned handle resolves at CALL time, so this is safe to call even when
% physics3d/ is not on the path; a 'shape3d' seed just needs physics3d/ on the
% path when it is actually integrated.
%
% INPUT   seedParams : struct with a .model field ('planar' | 'shape3d'). A struct
%                      without .model is treated as 'planar' (back-compatible).
% OUTPUT  rhs        : function handle f(t, x, seedParams) -> xdot.

    if isfield(seedParams, 'model')
        model = char(seedParams.model);
    else
        model = 'planar';
    end

    switch model
        case 'planar'
            rhs = @seed6DOFODE;
        case 'shape3d'
            rhs = @seed6DOFODE3D;
        otherwise
            error('seedRHS:unknownModel', ...
                  'Unknown seedParams.model "%s" (expected ''planar'' or ''shape3d'').', model);
    end
end
