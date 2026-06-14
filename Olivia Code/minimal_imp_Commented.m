% =========================================================================
%  2D falling / fluttering flat-plate dynamics (quasi-steady aero model)
%  ------------------------------------------------------------------------
%  This is the classic Andersen–Pesavento–Wang style model for a rigid 2D
%  plate moving through a fluid: a body-frame Newton–Euler set of equations
%  with anisotropic ADDED MASS, plus an empirical quasi-steady model for the
%  fluid lift, drag, center-of-pressure, and rotational torques.
%
%  All units are CGS: length in cm, mass in g, time in s
%  (so g = 980 cm/s^2 and water density = 1 g/cm^3).
% =========================================================================

clear;          % wipe workspace variables
clc;            % clear the command window
close all       % close any open figure windows

% --- Physical / geometric parameters (all CGS) ----------------------------
lcm = 0.1;      % offset of the center of mass (CoM) from the plate's
                %   geometric center, measured along the body x'-axis.
                %   Nonzero lcm is what makes the plate asymmetric and gives
                %   a buoyancy/gravity restoring couple.
l   = 2.54;     % plate chord length (full length, ~1 inch)
h   = 0.05;     % plate thickness (height)
g   = 980;      % gravitational acceleration (cm/s^2)
p   = 1;        % fluid (mass) density rho — 1 g/cm^3 = water
m   = 1;        % plate mass per unit depth (2D mass)
I   = .5;       % plate moment of inertia per unit depth (2D, about CoM)
mp  = p * l * h;% mass of fluid displaced by the plate cross-section.
                %   Used for the buoyancy correction: the EFFECTIVE weight
                %   driving the motion is (m - mp)*g, and mp also sets the
                %   buoyancy torque.

% --- Initial state -------------------------------------------------------
% State vector x = [ x ; y ; theta ; vxp ; vyp ; omega ]
%   x, y     : position of the GEOMETRIC center in the world frame
%   theta    : plate orientation angle (world)
%   vxp, vyp : translational velocity expressed in the BODY (primed) frame,
%              x' along the chord, y' normal to the plate
%   omega    : angular velocity (dtheta/dt)
x_0 = [0; 0; -pi/6; 5; 0; 0];   % start tilted -30 deg, moving at 5 cm/s along
                                %   the body x'-axis, no rotation

timerange = 5;                  % total simulated time (s)

% Tight integration tolerances — needed because this system is stiff-ish and
% sensitive (chaotic tumbling/fluttering), so loose tolerances drift badly.
option = odeset('RelTol',1e-9, 'AbsTol',1e-8);   % (defined but not passed to ode below)

% Build the ODE object (newer MATLAB ode/solve interface, R2023b+).
% The right-hand side is Model(...); extra params are bound via the anon fn.
F = ode(ODEFcn=@(t,x) Model(x,lcm,I,l,m,mp,p,g), InitialValue=x_0);

dt     = 0.001;                 % output sampling timestep (s)
Ntimes = timerange/dt;          % number of output samples
t      = linspace(0,timerange,Ntimes);  % requested output time grid

sol = solve(F,t);               % integrate the ODE on that time grid

% --- Extract and plot the trajectory -------------------------------------
t        = sol.Time;            % actual solution times
solution = sol.Solution;        % 6 x N state history
x = solution(1,:);              % geometric-center x(t)
y = solution(2,:);              % geometric-center y(t)
plot(x,y);                      % plot the flight path in the x-y plane
xlabel('x');
ylabel('y');
axis equal
grid on

% =========================================================================
function y = Model(initial,lcm,I,l,m,mp,p,g)
% Returns y = d/dt of the state vector:
%   y = [ dot_x, dot_y, dot_theta, dot_vxp, dot_vyp, dot_omega ]
% Note: x,y are the GEOMETRIC center; the primed (x',y') velocities and the
% dynamics are referenced to the CENTER OF MASS, offset by lcm along x'.

    y = zeros(size(initial));   % preallocate the derivative vector

    % --- Empirical aerodynamic-model coefficients ------------------------
    % These are fit/tuned constants for the quasi-steady flat-plate model.
    CL1  = 5.2;     % attached-flow lift slope (thin-airfoil-like; ~2*pi=6.28 for ideal)
    CL2  = 0.95;    % separated/stalled-flow lift amplitude (sin(2*alpha) term)
    CD0  = 0.1;     % parasitic (zero-incidence) drag coefficient
    CD1  = 5;       % attached-flow induced-drag coefficient (~sin^2 alpha)
    CD2  = 1.9;     % separated-flow (broadside) drag coefficient
    CR   = 1.1;     % rotational-lift / added-circulation coefficient (omega term)
    CCP0 = 0.3;     % center-of-pressure model: baseline offset (attached)
    CCP1 = 3.5;     % center-of-pressure model: curvature vs angle (attached)
    CCP2 = 0.2;     % center-of-pressure model: high-angle offset (separated)

    % --- Added (virtual) mass and inertia of the entrained fluid ---------
    m1 = 0;                       % added mass m_11: motion along the chord
                                  %   (edge-on) -> negligible for a thin plate
    m2 = pi * p * l^2 * 0.25;     % added mass m_22 = pi*rho*(l/2)^2: broadside
                                  %   motion entrains a fluid cylinder of dia. l
    Ia = pi * p * l^4 * (1 + 32*lcm^2/l^2)/128;
                                  % added moment of inertia about the CoM.
                                  %   Base term pi*rho*l^4/128 (rotational added
                                  %   inertia of the plate), plus a parallel-axis
                                  %   correction: note (32*lcm^2/l^2)*base = m2*lcm^2,
                                  %   i.e. it adds the added-mass term m2*lcm^2 for
                                  %   the CoM being offset from the geometric center.

    % --- Unpack the current state ----------------------------------------
    c     = cos(initial(3));      % cos(theta)
    s     = sin(initial(3));      % sin(theta)
    vxp   = initial(4);           % body-frame x'-velocity
    vyp   = initial(5);           % body-frame y'-velocity
    omega = initial(6);           % angular velocity

    % Velocity AT THE CENTER OF MASS. Because the CoM is offset by lcm along
    % x', rotation contributes -omega*lcm to the normal (y') velocity there.
    v_cg  = sqrt(vxp^2 + (vyp - omega*lcm)^2);   % speed of the CoM
    alpha = atan2(vyp - omega*lcm, vxp);         % angle of attack of the CoM
                                                 %   velocity relative to the chord

    % --- Smooth blending (weighting) functions ---------------------------
    % tanh switches smoothly between the ATTACHED model (small |alpha|) and the
    % SEPARATED model (large |alpha|) around a stall angle "as" with width "d".
    as = 14*pi/180;   % stall angle ~14 deg
    d  =  6*pi/180;   % blend width ~6 deg
    f1 = ( 1 - tanh((abs(alpha)-as)/d) )/2;      % ~1 attached, ~0 stalled (|alpha|<pi/2)
    f2 = ( 1 - tanh((pi-abs(alpha)-as)/d) )/2;   % same blend mirrored near alpha=+/-pi

    % --- Lift / normal-force coefficient (here labeled CT) ---------------
    % Blend of attached (CL1*sin(alpha)) and separated (CL2*sin(2*alpha)) lift.
    % Three branches cover the three angle ranges; the sign flips set the
    % correct lift direction front-vs-back of the plate.
    CT_1 = CL1*f1*sin(alpha)         + (1-f1)*CL2*sin(2*alpha);                 % alpha in (-90,90) deg
    CT_2 = -( CL1*f2*(pi-abs(alpha)) + (1-f2)*CL2*sin(2*(pi-abs(alpha))) );     % alpha in ( 90,180) deg
    CT_3 =    CL1*f2*(pi-abs(alpha)) + (1-f2)*CL2*sin(2*(pi-abs(alpha)));       % alpha in (-180,-90) deg
    % NOTE: in CT_2/CT_3 the first term uses (pi-abs(alpha)) directly, whereas
    % CT_1 uses sin(alpha). This looks like a small-angle approximation of
    % sin(pi-|alpha|) ~ (pi-|alpha|) near alpha = +/-pi; worth confirming it is
    % intentional rather than a missing sin().

    % --- Drag coefficient ------------------------------------------------
    % Attached: CD0 + CD1*sin^2(alpha); separated: CD2*sin^2(alpha). Blended.
    CD_1 = f1*(CD0 + CD1*(sin(alpha))^2)        + (1-f1)*CD2*(sin(alpha))^2;            % |alpha|<90
    CD_2 = f2*(CD0 + CD1*(pi-abs(alpha))^2)     + (1-f2)*CD2*(sin(pi-alpha))^2;         % |alpha|>90
    % NOTE (likely inconsistency): CD_2 uses (pi-abs(alpha))^2 in the first term
    % (vs sin^2 in CD_1) and sin(pi-alpha) (no abs) in the second term. For
    % consistency these probably should be sin(pi-abs(alpha))^2; flag to verify.

    % --- Center of pressure (distance from geometric center) -------------
    r  = abs(alpha);    % |alpha|
    rr = pi - r;        % supplement of |alpha| (used for the high-angle branch)
    % Attached branch: quadratic-in-angle CoP; separated branch: linear ramp.
    lcp_1 =  f1*(CCP0 - CCP1*r^2)  + (1-f1)*(CCP2 - (CCP2/(pi/2))*r);   % low-angle CoP fraction
    lcp_1 =  lcp_1*l;                                                   % -> physical length
    lcp_2 = -( f2*(CCP0 - CCP1*rr^2) + (1-f2)*(CCP2 - (CCP2/(pi/2))*rr) );% high-angle CoP fraction (sign-flipped)
    lcp_2 =  lcp_2*l;                                                   % -> physical length

    % --- Select the correct branch for the current angle of attack -------
    if r<=pi/2              % |alpha| <= 90 deg  : attached/forward regime
        CT_alpha = CT_1;
        CD_alpha = CD_1;
        l_cp     = lcp_1;
    elseif alpha>pi/2       % alpha in (90,180)  : reversed flow, front
        CT_alpha = CT_2;
        CD_alpha = CD_2;    % drag is symmetric -> reuse CD_2 / lcp_2
        l_cp     = lcp_2;
    else                    % alpha in (-180,-90): reversed flow, back
        CT_alpha = CT_3;    % only the lift sign differs from the elseif branch
        CD_alpha = CD_2;
        l_cp     = lcp_2;
    end

    l_tau = l_cp - lcm;     % moment arm of the aero force about the CoM
                            %   (CoP location minus CoM offset)
    lcr   = 0;              % reference point for the rotational-lift torque
                            %   (here taken at the geometric center)

    % --- Lift force (perpendicular to the relative velocity) -------------
    % L_0 bundles the lift coefficient * one power of speed, PLUS a rotational
    % (added-circulation, "-l*CR*omega") contribution. Multiplying by the
    % velocity components below supplies the 2nd speed power and the direction.
    L_0 = 0.5*p*l*(CT_alpha*v_cg - l*CR*omega);
    Lx  =  L_0*(vyp - omega*lcm);   % x'-component of lift (perp. to velocity)
    Ly  = -L_0*vxp;                 % y'-component of lift (perp. to velocity)

    % --- Drag force (opposes the relative velocity) ----------------------
    D_0  = -0.5*p*l*CD_alpha*v_cg;  % drag magnitude factor (negative = opposing)
    D_xp = D_0*vxp;                 % x'-component of drag (along -velocity)
    D_yp = D_0*(vyp - omega*lcm);   % y'-component of drag (along -velocity)

    % --- Torques about the center of mass --------------------------------
    % Torque from the translational aero force acting at the CoP (force*arm):
    Tt  = -0.5*p*l*v_cg*(CT_alpha*vxp + CD_alpha*(vyp - omega*lcm))*l_tau;
    % Torque from the rotational lift (the CR / omega term), arm (lcm-lcr):
    Tcr = -0.5*p*l^2*CR*omega*vxp*(lcm - lcr);
    % Rotational (spin) DAMPING torque ~ omega*|omega| (quadratic in rate),
    % from integrating the pressure distribution along the plate about the
    % offset CoM -> the ((2lcm/l +/- 1)^4) terms are the two plate-end limits.
    Tr  = -1/128*p*l^4*CD2*omega*abs(omega)*((2*lcm/l + 1)^4 + (2*lcm/l - 1)^4);
    % Buoyancy/gravity restoring couple: buoyant force acts at the geometric
    % center, weight at the CoM (offset lcm); cos(theta) projects the lever arm.
    Tb  = -mp*g*lcm*c;

    % --- Assemble equations of motion (body frame, with added mass) ------
    total_torque = Tt + Tr + Tb + Tcr;   % net torque about the CoM

    % Body-frame force balances. The omega*v terms are the rotating-frame
    % (Coriolis/centripetal) coupling, made anisotropic by the added masses.
    % The -m2*omega^2*lcm and +m2*(domega/dt)*lcm terms come from the CoM being
    % offset from the geometric center. The last term is the effective weight
    % (m-mp)*g resolved into the body axes via s = sin(theta), c = cos(theta).
    F_xp =  (m+m2)*omega*vyp - m2*omega^2*lcm        + Lx + D_xp - (m-mp)*g*s;
    F_yp = -(m+m1)*omega*vxp + m2*total_torque/(I+Ia)*lcm + Ly + D_yp - (m-mp)*g*c;
    % NOTE: F_yp reuses total_torque/(I+Ia) as a stand-in for the angular
    % acceleration (domega/dt) — i.e. the otherwise coupled (vyp,omega)
    % equations are solved explicitly/sequentially rather than as one linear
    % system. Fine if the coupling is weak, but worth being aware of.

    % --- Pack the state derivatives --------------------------------------
    y(1) = vxp*c - vyp*s;          % dx/dt  : body -> world rotation of velocity
    y(2) = vxp*s + vyp*c;          % dy/dt
    y(3) = omega;                  % dtheta/dt
    y(4) = F_xp / (m+m1);          % dvxp/dt  (divide by effective x'-mass)
    y(5) = F_yp / (m+m2);          % dvyp/dt  (divide by effective y'-mass)
    y(6) = total_torque / (I+Ia);  % domega/dt (divide by effective inertia)
end