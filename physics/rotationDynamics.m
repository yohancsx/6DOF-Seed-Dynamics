function alpha_body = rotationDynamics(tau_body, omega_body, mp)
% ROTATIONDYNAMICS  Angular acceleration in the body frame (modified Euler eqn).
%
% Solves Euler's equation about the CoM with added inertia and a possibly
% time-varying inertia tensor:
%   (I_G + A_rot) * alpha = tau - omega x ((I_G + A_rot)*omega) - I_G_dot*omega
%
%   - omega x (I_eff*omega) : gyroscopic term (identically zero in the scalar 2D
%                             case; real in 3D).
%   - I_G_dot*omega         : correction when the inertia changes in time (e.g. a
%                             drifting CoM). Zero for a rigid, CoM-fixed seed.
% The added-inertia rate (A_rot_dot) is neglected (small, and ~0 for a fixed CoM).
%
% INPUTS
%   tau_body   : 3x1 net torque about the CoM, body frame (N·m).
%   omega_body : 3x1 angular velocity, body frame (rad/s).
%   mp         : mass-properties struct (uses mp.I_G, mp.A_rot, mp.I_G_dot).
%
% OUTPUT
%   alpha_body : 3x1 angular acceleration, body frame (rad/s^2).

    I_eff = mp.I_G + mp.A_rot;                  % effective inertia (rigid + added)

    alpha_body = I_eff \ ( tau_body ...
                           - cross(omega_body, I_eff * omega_body) ...   % gyroscopic
                           - mp.I_G_dot * omega_body );                  % varying inertia
end
