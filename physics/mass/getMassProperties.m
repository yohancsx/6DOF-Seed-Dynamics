function mp = getMassProperties(t, seedParams)
% GETMASSPROPERTIES  Mass properties at the current simulation time.
%
% Interpolates the precomputed time series from setupSeedShapeAndMass
% (massParams.com_t, I_G_t, I_G_dot_t) at time t. The sample times
% (massParams.tSamples) are assumed to be the same time axis the simulation
% runs on. Behaviour at the edges:
%   - t beyond the last sample  -> clamped to the last value (CoM, I_G, I_G_dot)
%   - t before the first sample -> clamped to the first value (defensive)
%   - a single sample           -> treated as constant (no interpolation)
%
% INPUT
%   t          : current simulation time (s), scalar.
%   seedParams : struct from setupSeedShapeAndMass, using sub-struct
%                seedParams.massParams with fields:
%                  .tSamples  (1xN or Nx1) sample times (s)
%                  .com_t     (3xN)        CoM in body-datum coords (m)
%                  .I_G_t     (3x3xN)      inertia about CoM (kg·m^2)
%                  .I_G_dot_t (3x3xN)      d/dt of inertia (kg·m^2/s)
%                  .M_total   (scalar)     total mass (kg)
%
% OUTPUT (struct mp)
%   mp.M       : total mass (kg), scalar.
%   mp.c       : 3x1 CoM in body-datum coords at time t (m).
%   mp.I_G     : 3x3 inertia about the CoM at time t (kg·m^2).
%   mp.I_G_dot : 3x3 inertia rate at time t (kg·m^2/s).
%   mp.A_trans : 3x3 translational added-mass block, body frame (kg).
%   mp.A_rot   : 3x3 rotational added-inertia block about the CoM (kg·m^2).
%
% Added mass is evaluated at the CURRENT CoM (mp.c) via getAddedMass, so a
% moving CoM is tracked. Buoyancy geometry (Vdisp, centroid) is intentionally
% omitted -- buoyancy is negligible for seeds.

% -------------------------------------------------------------------------
% 0. UNPACK
% -------------------------------------------------------------------------
massParams = seedParams.massParams;
tSamples   = massParams.tSamples(:);   % force column for indexing/interp
numT       = numel(tSamples);

% -------------------------------------------------------------------------
% 1. TOTAL MASS  (constant)
% -------------------------------------------------------------------------
mp.M = massParams.M_total;

% -------------------------------------------------------------------------
% 2. CoM AND INERTIA AT TIME t
% -------------------------------------------------------------------------
if numT == 1
    % Constant case: a single time sample -> no interpolation.
    mp.c       = massParams.com_t(:, 1);
    mp.I_G     = massParams.I_G_t(:, :, 1);
    mp.I_G_dot = massParams.I_G_dot_t(:, :, 1);
else
    % Clamp the query time to the sampled span before interpolating: t past
    % the end returns the last value, t before the start the first, no NaNs.
    tq = min(max(t, tSamples(1)), tSamples(end));

    % CoM: com_t is 3xN -> interp each component over time.
    mp.c = interp1(tSamples, massParams.com_t.', tq).';   % 3x1

    % Inertia tensors: flatten 3x3xN to Nx9, interpolate, reshape to 3x3.
    IG_flat    = reshape(massParams.I_G_t,     9, numT).';   % Nx9
    IGdot_flat = reshape(massParams.I_G_dot_t, 9, numT).';   % Nx9

    mp.I_G     = reshape( interp1(tSamples, IG_flat,    tq), 3, 3 );
    mp.I_G_dot = reshape( interp1(tSamples, IGdot_flat, tq), 3, 3 );
end

% -------------------------------------------------------------------------
% 3. ADDED MASS (evaluated at the current CoM)
%    A_trans is CoM-independent; A_rot is referenced to the CoM, so it is
%    re-evaluated at mp.c each call and tracks a moving CoM.
% -------------------------------------------------------------------------
[mp.A_trans, mp.A_rot] = getAddedMass(seedParams, mp.c);

end