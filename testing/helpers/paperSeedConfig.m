function [cfg, bsp, toNutPos, info] = paperSeedConfig()
% PAPERSEEDCONFIG  The Hou et al. (2025) samara-inspired plate, in our units.
%
% One place that builds THEIR seed, so the anchor-case script and the paper mode
% grid cannot drift apart. Source: Hou, Zhang, Li, Jia & Huang, "Aerodynamic
% significance of mass distribution on diverse samara descent behaviors",
% Commun. Eng. 4, 129 (2025), Methods > "Samara-inspired framework".
%
% =========================================================================
% THEIR RIG, AND HOW EACH NUMBER MAPS ONTO OURS
% =========================================================================
%   Plate      60 mm x 20 mm kraft paper, aspect ratio 3.0, mass 92.8 mg.
%              Their LONGER axis is our SPAN (body z); their shorter axis is our
%              CHORD (body x). (setupSeedShapeAndMass reads the polyshape's x as
%              body z and its y as body x, so the polyshape is built that way.)
%   Weights    a heavier 102.9 mg on the longer axis, a lighter 51.4 mg on the
%              shorter one. WE LUMP THEM into a single 154.3 mg nut, which
%              reproduces their CoM and their TOTAL mass (247.1 mg) exactly, and
%              therefore their average surface density and descent-velocity
%              scaling. It does NOT reproduce their inertia: their light weight
%              sits up to 10 mm off the long axis and contributes ~4.3e-9 kg m^2
%              about it, where our single equivalent nut contributes ~8.2e-10 --
%              5x less, about the very axis spiral tumbling rotates around. That
%              is the phase-5 Tier 2 item; treat a wrong ST tumble frequency as
%              this, not as physics.
%   Thickness  NOT STATED in the paper. Derived: their plate's areal density is
%              92.8 mg / (60 x 20 mm) = 77.3 g/m^2, and kraft paper of that
%              grammage runs about 0.10 mm thick (bulk density ~770 kg/m^3), so
%              t = 1.0e-4 m. It matters for exactly one thing in this model, the
%              edge-drag frontal area t*c (setupSeedShapeAndMass reads thickness
%              for buoyancy volume but never uses it). With their thin plate that
%              area is ~15x smaller than our own 2 mm seed's, so a span-first fall
%              is resisted ~4x less -- if a case collapses edge-on it will run away
%              to tens of m/s. That is a DIAGNOSTIC, not a bug: it says tip form
%              drag alone is too weak a spanwise resistance at this thickness.
%   Air        Experiments at 15 +/- 1 C; rho = 1.225 kg/m^3 is right for that.
%   Release    Their clamp releases at cone angle theta = 0 (long axis horizontal)
%              and self-rotation psi = pi/2 (plate face vertical, edge-on to
%              gravity), from rest. Returned as cfg.releaseQuat / cfg.releaseOmega.
%              NOTE this differs from our own suites, which release at a 30 deg
%              tilt about the span -- a different basin of attraction.
%   Drop       4 m enclosed platform; cfg.dropHeight records it so a run can be
%              read at their observation window as well as at settling time.
%
% =========================================================================
% THE PARAMETER SPACE (their Fig. 2a axes)
% =========================================================================
%   x_c = m_h * x'/m_total   along the LONGER axis,  a = L/2  ->  x_c/a <= m_h/m_total = 0.416
%   y_c = m_l * y'/m_total   along the SHORTER axis, b = W/2  ->  y_c/b <= m_l/m_total = 0.208
% Those caps are structural, and they are exactly where their plotted data stops,
% which is what confirms the reading of the axes. For us, with a single nut of
% mass fraction mu = m_nut/M_total,
%       x_c/a = mu * z_nut/(L/2)        y_c/b = mu * x_nut/(W/2)
% inverted by toNutPos below. With mu = 0.625 the whole window is reachable with
% the nut ON the plate (z_nut <= 20 mm of a 30 mm half-span).
%
% OUTPUTS
%   cfg      : config struct for buildSeedParams / runSingleMode (shape3d).
%   bsp      : base-seed-params for their plate (nut position set by the caller,
%              or by runSingleMode from cfg.nutMass).
%   toNutPos : @(xc_over_a, yc_over_b) -> 3x1 body nut position [x;0;z] (m).
%   info     : struct of derived quantities -- .L .W .mu .plateMass .nutMass
%              .totalMass .sigma .Vscale (= sqrt(sigma*g/rho)) .arealDensity
%              .xcMax .ycMax (the structural caps above).

    % --- Their rig ---------------------------------------------------------
    L         = 0.060;        % plate length  -> our SPAN  (body z), m
    W         = 0.020;        % plate width   -> our CHORD (body x), m
    plateMass = 92.8e-6;      % kraft plate, kg
    m_h       = 102.9e-6;     % heavier weight, on the longer axis, kg
    m_l       = 51.4e-6;      % lighter weight, on the shorter axis, kg

    arealDensity = plateMass / (L * W);       % 0.0773 kg/m^2 (= 77.3 g/m^2)
    thickness    = 1.0e-4;                    % derived, see header
    nutMass      = m_h + m_l;                 % lumped, see header
    totalMass    = plateMass + nutMass;
    mu           = nutMass / totalMass;
    sigma        = totalMass / (L * W);       % their "average surface density"

    % --- Config -------------------------------------------------------------
    cfg.spanLength  = L;
    cfg.chordLength = W;
    cfg.thickness   = thickness;
    cfg.bulkDensity = arealDensity / thickness;   % so bulkDensity*thickness = areal
    cfg.numStrips   = 12;
    cfg.tSamples    = 0;
    cfg.nutMass     = nutMass;
    cfg.rhoFluid    = 1.225;
    cfg.g           = 9.81;
    cfg.shapeModel  = 'shape3d';                  % flat: they excluded camber and twist

    cfg.tspan      = [0 12];                      % long enough to settle; their 4 m
    cfg.dropHeight = 4.0;                         % fall is ~3 s at these speeds
    cfg.odeRelTol  = 1e-6;
    cfg.odeAbsTol  = 1e-8;
    cfg.metricOpts.windowStartFrac = 0.5;
    cfg.metricOpts.convergeTol     = 0.20;
    cfg.modeThresholds = defaultModeThresholds();

    % Release: theta = 0, psi = pi/2 -- long axis horizontal, face vertical.
    cfg.releaseQuat  = axisAngleToQuat([0; 0; 1], pi/2);
    cfg.releaseOmega = [0; 0; 0];

    % Physics switches deliberately left UNSET, so this config always tracks the
    % shape3d builder's defaults rather than pinning a stale copy of them. As of
    % the Tier 1 screen those are edge drag ON, both halves of Kirchhoff's pair ON
    % (added-mass rate and Munk moment), LEV OFF. Override through cfg, not here.

    % --- Base seed params ---------------------------------------------------
    hL = L/2;   hW = W/2;
    bsp.seedShape     = polyshape([-hL, hL, hL, -hL], [-hW, -hW, hW, hW]);
    bsp.seedDensity   = arealDensity;
    bsp.seedThickness = thickness;
    bsp.numStrips     = cfg.numStrips;
    bsp.tSamples      = 0;
    bsp.nutPos_t      = [0; 0; 0];
    bsp.nutMass_t     = nutMass;

    % --- (x_c/a, y_c/b) -> nut position ------------------------------------
    toNutPos = @(xa, yb) [ yb * hW / mu ; 0 ; xa * hL / mu ];

    % --- Derived quantities worth having at hand ---------------------------
    info = struct('L', L, 'W', W, 'plateMass', plateMass, 'nutMass', nutMass, ...
                  'm_h', m_h, 'm_l', m_l, 'totalMass', totalMass, 'mu', mu, ...
                  'arealDensity', arealDensity, 'thickness', thickness, ...
                  'sigma', sigma, 'Vscale', sqrt(sigma * cfg.g / cfg.rhoFluid), ...
                  'xcMax', m_h / totalMass, 'ycMax', m_l / totalMass);
end
