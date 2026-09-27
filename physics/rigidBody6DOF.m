function [xdot, core] = rigidBody6DOF(x, mp, F_aero_body, tau_body, gMag)
% RIGIDBODY6DOF  Shape-agnostic 6-DOF rigid-body state derivative.
%
% The part of the seed equations of motion that does NOT depend on the seed's
% shape. Given the current state, the mass properties (including added mass), and
% the NET aerodynamic force/torque already summed in the BODY frame, it rotates
% the force to inertial, adds gravity, solves the translational dynamics (with
% translational added mass) and the modified-Euler rotational dynamics (with I_G,
% the I_G_dot transport term, and rotational added mass), integrates the
% quaternion, and assembles the 13-state derivative.
%
% Both the planar force model (seed6DOFODE) and the future 3D-shape force model
% feed this the SAME way: each computes its own (F_aero_body, tau_body) and mass
% properties, then calls this core. Keeping the core here means the dynamics that
% must never differ between the planar and 3D models live in exactly one place.
%
% INPUTS
%   x           : 13x1 state [ r(3); q(4); v(3); omega(3) ].
%                   q     body->world quaternion, scalar-first (normalised here)
%                   v     CoM velocity, INERTIAL frame (m/s)
%                   omega angular velocity, BODY frame (rad/s)
%                 (r is not used for the derivative.)
%   mp          : mass-properties struct from getMassProperties, with fields
%                 .M .c .I_G .I_G_dot .A_trans .A_rot (as consumed by
%                 translationDynamics / rotationDynamics).
%   F_aero_body : 3x1 net aerodynamic force,  BODY frame (N).
%   tau_body    : 3x1 net aerodynamic torque, BODY frame (N*m), about the CoM.
%   gMag        : gravity magnitude (m/s^2); gravity acts along world -Y (Y-up).
%
% OUTPUTS
%   xdot : 13x1 state derivative [ v; dq; a_inertial; alpha_body ].
%   core : (optional) shared intermediates for post-processing:
%          .F_aero_inertial, .F_total_inertial (aero + gravity, as before),
%          .F_addedMassRate (the Adot*v term, zero when switched off),
%          .tau_addedMass   (the v x (A*v) Munk moment, zero when switched off),
%          .a_inertial, .alpha_body.
%
% ADDED-MASS RATE (mp.enableAddedMassRate, default FALSE)
% The seed drags a little fluid along with it, and how much -- and in which
% direction -- depends on which way the plate faces. The translational added mass
% is fixed in the BODY frame (A_b) but the body rotates, so in the inertial frame
% where translation is solved it is A_i = R*A_b*R'. The fluid's momentum is
% (M*I + A_i)*v, and Newton applies to its rate of change:
%     (M*I + A_i)*a = F - Adot_i*v
% With A_b constant, Rdot = [w x]*R (w = inertial angular velocity) gives
%     Adot_i*v = w x (A_i*v) - A_i*(w x v)
% i.e. spinning changes the entrained fluid's momentum even at constant speed.
% translationDynamics solves (M*I + A_i)*a = F; this term is folded into the F it
% receives, so its signature is unchanged. The formula was verified against a
% finite difference of R*A_b*R' to 6.2e-9, and measured at 31-43% of the net force
% in spinning cases -- not the "tiny" term earlier docstrings assumed.
%
% Provenance: this is Kirchhoff's equation for a body moving in a fluid with an
% anisotropic added mass (Lamb, Hydrodynamics, ch. VI), written in the inertial
% frame. The APW 2D reference model carries the SAME physics in its body-frame
% form -- the (m+m2)*omega*vyp and -(m+m1)*omega*vxp terms of minimal_imp
% (Olivia Code/minimal_imp_Commented.m, lines 208-209). Written in the body frame
% those cross terms appear explicitly; written in the inertial frame they become
% this Adot*v. testNewPhysicsTerms check A5 shows the two forms agree to machine
% precision with the term ON and disagree with it OFF.
%
% It is OFF by default in BOTH builders (planar, to stay byte-identical; shape3d,
% by decision while its support is reviewed). Set enableAddedMassRate = true to
% include it. (The rotational counterpart Adot_rot*omega is always neglected:
% measured at 0.03-0.1% of the applied torque.)
%
% ADDED-MASS MOMENT / MUNK MOMENT (mp.enableAddedMassMoment, default FALSE)
% The other half of Kirchhoff's equations, and the ROTATIONAL partner of the term
% above. Kirchhoff's pair for a body with added mass A in an inviscid fluid (Lamb,
% Hydrodynamics ch. VI), body frame, with p = (M*I + A)*v and L = (I_G + A_rot)*w:
%     dp/dt + w x p = F
%     dL/dt + w x L + v x p = tau
% The v x p term is the one this switch adds. Since v x (M*v) = 0, it is
%     v x (A*v)
% and it moves to the right-hand side as -v x (A*v), exactly as Adot_i*v does for
% the force. It is NOT a small correction: it is the moment that turns a body
% broadside-on to its own motion, i.e. what makes a plate falling EDGE-ON
% unstable. With A anisotropic (A_normal >> A_chord, A_span for a plate) it is
% quadratic in velocity and needs no coefficient -- it comes from the added-mass
% tensor the model already carries.
%
% Provenance: the APW reference model carries it explicitly. Andersen, Pesavento &
% Wang (2005), J. Fluid Mech. 541:65-90, eq. (6.3):
%     (I + I_a) thetaddot = (m11 - m22) vx' vy'  +  l_tau rho_f Gamma |v|  -  tau_visc
% Our strip forces supply their second term (circulatory torque at l_cp) and their
% third (stripSpinDamping), but the FIRST term had no counterpart here at all. In
% 2D this implementation reduces to it exactly: with A = diag(m11, m22, .) and
% v = (vx', vy', 0),  -[v x (A v)]_z = (m11 - m22) vx' vy'  (testNewPhysicsTerms
% check M2). Their m11 = (pi/4) rho_f h^2 and m22 = (pi/4) rho_f l^2 are the same
% inviscid coefficients getAddedMass builds.
%
% OFF by default in both builders, so nothing changes until it is switched on.

    % --- Parse state (normalise q defensively; state may have drifted) ------
    q     = x(4:7);   q = q / norm(q);
    v     = x(8:10);          % CoM velocity, inertial
    omega = x(11:13);         % angular velocity, body frame
    R     = quatToRotm(q);    % body->world

    % --- Forces: aero (body->inertial) + gravity ---------------------------
    F_aero_inertial  = R * F_aero_body;
    F_grav           = mp.M * [0; -gMag; 0];        % world -Y; buoyancy neglected
    F_total_inertial = F_aero_inertial + F_grav;

    % --- Added-mass rate:  F_solve = F_total - Adot_i*v  (see header) -------
    F_addedMassRate = [0; 0; 0];
    if isfield(mp, 'enableAddedMassRate') && mp.enableAddedMassRate
        A_i     = R * mp.A_trans * R.';              % added mass, inertial frame
        omega_i = R * omega;                         % angular velocity, inertial
        F_addedMassRate = cross(omega_i, A_i*v) - A_i*cross(omega_i, v);   % = Adot_i*v
    end
    F_solve = F_total_inertial - F_addedMassRate;

    % --- Added-mass (Munk) moment:  tau_solve = tau - v x (A*v)  (see header) -
    % Formed in the BODY frame, where A_trans lives and where the torque is
    % applied, so no rotation of the tensor is needed.
    tau_addedMass = [0; 0; 0];
    if isfield(mp, 'enableAddedMassMoment') && mp.enableAddedMassMoment
        v_body        = R.' * v;
        tau_addedMass = cross(v_body, mp.A_trans * v_body);      % = v x (A v)
    end
    tau_solve = tau_body - tau_addedMass;

    % --- Accelerations ------------------------------------------------------
    a_inertial = translationDynamics(F_solve, R, mp);            % with added mass
    alpha_body = rotationDynamics(tau_solve, omega, mp);         % modified Euler

    % --- Orientation kinematics --------------------------------------------
    dq = quatKinematics(q, omega);

    % --- Assemble (dr/dt = v, since inertial velocity is a state) -----------
    xdot = [ v ; dq ; a_inertial ; alpha_body ];

    if nargout > 1
        core.F_aero_inertial  = F_aero_inertial;
        core.F_total_inertial = F_total_inertial;
        core.F_addedMassRate  = F_addedMassRate;
        core.tau_addedMass    = tau_addedMass;
        core.a_inertial       = a_inertial;
        core.alpha_body       = alpha_body;
    end
end
