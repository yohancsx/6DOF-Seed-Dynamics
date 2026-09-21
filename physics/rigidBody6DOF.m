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
% in spinning cases -- not the "tiny" term earlier docstrings assumed. It is OFF by
% default so the frozen planar model stays byte-identical; setupSeedShape3D turns
% it on for the shape3d model. (The rotational counterpart Adot_rot*omega is still
% neglected: measured at 0.03-0.1% of the applied torque.)

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

    % --- Accelerations ------------------------------------------------------
    a_inertial = translationDynamics(F_solve, R, mp);            % with added mass
    alpha_body = rotationDynamics(tau_body, omega, mp);          % modified Euler

    % --- Orientation kinematics --------------------------------------------
    dq = quatKinematics(q, omega);

    % --- Assemble (dr/dt = v, since inertial velocity is a state) -----------
    xdot = [ v ; dq ; a_inertial ; alpha_body ];

    if nargout > 1
        core.F_aero_inertial  = F_aero_inertial;
        core.F_total_inertial = F_total_inertial;
        core.F_addedMassRate  = F_addedMassRate;
        core.a_inertial       = a_inertial;
        core.alpha_body       = alpha_body;
    end
end
