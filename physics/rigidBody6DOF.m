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
%          .F_aero_inertial, .F_total_inertial, .a_inertial, .alpha_body.

    % --- Parse state (normalise q defensively; state may have drifted) ------
    q     = x(4:7);   q = q / norm(q);
    v     = x(8:10);          % CoM velocity, inertial
    omega = x(11:13);         % angular velocity, body frame
    R     = quatToRotm(q);    % body->world

    % --- Forces: aero (body->inertial) + gravity ---------------------------
    F_aero_inertial  = R * F_aero_body;
    F_grav           = mp.M * [0; -gMag; 0];        % world -Y; buoyancy neglected
    F_total_inertial = F_aero_inertial + F_grav;

    % --- Accelerations ------------------------------------------------------
    a_inertial = translationDynamics(F_total_inertial, R, mp);   % with added mass
    alpha_body = rotationDynamics(tau_body, omega, mp);          % modified Euler

    % --- Orientation kinematics --------------------------------------------
    dq = quatKinematics(q, omega);

    % --- Assemble (dr/dt = v, since inertial velocity is a state) -----------
    xdot = [ v ; dq ; a_inertial ; alpha_body ];

    if nargout > 1
        core.F_aero_inertial  = F_aero_inertial;
        core.F_total_inertial = F_total_inertial;
        core.a_inertial       = a_inertial;
        core.alpha_body       = alpha_body;
    end
end
