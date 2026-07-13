function [A_trans, A_rot] = getAddedMass(seedParams, comPos_body)
%TODO: CHECK THIS
% GETADDEDMASS  Body-frame added-mass blocks for the flat-plate seed.
%
% Strip-theory generalization of the 2D Andersen-Pesavento-Wang added mass
% (minimal_imp_Commented.m lines 92-102), which used:
%     m1 = 0                         % along-chord (edge-on) added mass
%     m2 = pi*rho*(l/2)^2            % broadside (normal) added mass
%     Ia = pi*rho*l^4/128 * (1 + 32*lcm^2/l^2)
%        = pi*rho*l^4/128  +  m2*lcm^2   % local term + parallel-axis to the CoM
% Here each spanwise strip contributes its own 2D added mass, and these are
% summed (with the proper 3D lever arms) into a translational block A_trans and
% a rotational block A_rot about the CoM. Buoyancy / displaced volume is NOT
% handled here (negligible for seeds, ignored by request).
%
% REASONING (why these specific terms):
%   * A thin flat plate moving EDGE-ON or IN ITS OWN PLANE entrains essentially
%     no fluid, so the chordwise (body x) and spanwise (body z) translational
%     added masses are ~0 -- exactly the 2D m1 = 0, extended to the new span
%     direction.
%   * Moving BROADSIDE (along the plate normal, body y) pushes a roughly
%     cylindrical slug of fluid of diameter = chord. Per unit span that 2D added
%     mass is m2 = pi*rho*(chord/2)^2; multiplying by the strip width dz_i gives
%     the strip's normal added mass. Summing over strips gives the total.
%   * For ROTATION, only motions that produce a NORMAL (body y) velocity excite
%     the (dominant) normal added mass m2; in-plane motions excite only the ~0
%     in-plane added masses. Working out omega x r for each body axis:
%       - spin about body z (spanwise): a chordwise offset x_s gives v_y = w*x_s
%         -> excites m2 with lever x_s, AND spins the cross-section in its own
%         plane -> excites the local 2D rotational added inertia ia. This is the
%         2D case: I_zz = ia + m2*x_s^2  (matches Ia = base + m2*lcm^2).
%       - spin about body x (chordwise): a spanwise offset z_s gives v_y = w*z_s
%         -> excites m2 with lever z_s:  I_xx = m2*z_s^2.  [NEW 3D]
%       - spin about body y (normal): sweeps the plate in-plane only -> excites
%         just the ~0 in-plane added masses:  I_yy ~ 0.  [NEW 3D]
%     The x-z cross term (-m2*x_s*z_s) is the product-of-inertia of the normal
%     added mass; it is the natural extension of the 2D parallel-axis term and
%     vanishes for a symmetric, CoM-centered seed.
%
% INPUTS
%   seedParams  : struct from setupSeedShapeAndMass. Uses:
%                   .strips.chord    (1xM) local chord per strip (m)
%                   .strips.dz       (1xM) strip spanwise width (m)
%                   .strips.xgc_body (1xM) strip chordwise geometric centre, body x (m)
%                   .strips.zgc_body (1xM) strip spanwise  geometric centre, body z (m)
%                 and the FLUID density seedParams.rhoFluid (kg/m^3). For a seed
%                 in air use ~1.225; added mass scales linearly with it.
%   comPos_body : (optional) 3x1 CoM in body-datum coords (m), the reference for
%                 A_rot. Defaults to the first column of massParams.com_t. For a
%                 moving CoM, pass the current value so A_rot tracks it.
%
% OUTPUTS
%   A_trans : 3x3 translational added-mass block, body frame (kg).
%             diag([0, sum(m2), 0]) -- only the normal (y) direction is nonzero.
%   A_rot   : 3x3 rotational added-inertia block about the CoM, body frame
%             (kg·m^2). Symmetric, with zero entries on the body-y row/column
%             except where coupling is physically absent.

% -------------------------------------------------------------------------
% 0. INPUTS / DEFAULTS
% -------------------------------------------------------------------------
if nargin < 2 || isempty(comPos_body)
    comPos_body = seedParams.massParams.com_t(:, 1);   % constant-CoM fallback
end

rhoFluid = seedParams.rhoFluid;        % fluid (air) density (kg/m^3)

chord = seedParams.strips.chord(:);    % Mx1 local chord (m)
dz    = seedParams.strips.dz(:);       % Mx1 strip width (m)
xGeo  = seedParams.strips.xgc_body(:); % Mx1 chordwise geometric centre, body x (m)
zGeo  = seedParams.strips.zgc_body(:); % Mx1 spanwise  geometric centre, body z (m)

% -------------------------------------------------------------------------
% 1. PER-STRIP 2D ADDED MASS / ADDED INERTIA
% -------------------------------------------------------------------------
% Normal (broadside) added mass of each strip: 2D m2 per unit span * width.
m2_strip = pi * rhoFluid .* (chord/2).^2 .* dz;        % Mx1 (kg)   [2D line 94]

% Local rotational added inertia of each strip cross-section about its own
% centroid (in-plane rotation about body z): 2D base term * width.
ia_strip = pi * rhoFluid .* chord.^4 .* dz / 128;      % Mx1 (kg·m^2) [2D line 96 base]

% -------------------------------------------------------------------------
% 2. TRANSLATIONAL ADDED-MASS BLOCK
%    Only the normal (body y) direction entrains fluid; chord (x) and span (z)
%    are edge-on/in-plane -> 0 (the 2D m1 = 0, extended to span).
% -------------------------------------------------------------------------
M_normal = sum(m2_strip);              % total broadside added mass (kg)
A_trans  = diag([0, M_normal, 0]);

% -------------------------------------------------------------------------
% 3. ROTATIONAL ADDED-INERTIA BLOCK (about the CoM)
%    Lever arms of each strip centroid relative to the CoM:
% -------------------------------------------------------------------------
x_s = xGeo - comPos_body(1);           % chordwise offset from CoM (body x)
z_s = zGeo - comPos_body(3);           % spanwise  offset from CoM (body z)

% Diagonal terms:
I_xx = sum( m2_strip .* z_s.^2 );                 % spin about chord axis  [NEW 3D]
I_zz = sum( ia_strip + m2_strip .* x_s.^2 );      % spin about span axis   [2D Ia]
I_yy = 0;                                         % spin about normal axis ~0 [NEW 3D]

% Off-diagonal x-z coupling (product of inertia of the normal added mass):
I_xz = -sum( m2_strip .* x_s .* z_s );            % vanishes if CoM-centered/symmetric

A_rot = [ I_xx,   0,    I_xz ;
          0,      I_yy, 0    ;
          I_xz,   0,    I_zz ];

end