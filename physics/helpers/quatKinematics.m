function dq = quatKinematics(q, omegaBody)
% QUATKINEMATICS  Time-derivative of the orientation quaternion.
%
% Advances a body->world quaternion given the BODY-frame angular velocity:
%       dq/dt = 1/2 * q (x) [0; omegaBody]
% A small normalisation-feedback term is added so the quaternion stays close to
% unit length under numerical integration (ode45 would otherwise let |q| drift).
%
% INPUTS
%   q         : 4x1 orientation quaternion [q0; q1; q2; q3], body->world.
%   omegaBody : 3x1 angular velocity in the BODY frame (rad/s).
%
% OUTPUT
%   dq : 4x1 quaternion time-derivative.

    % Pure quaternion built from the body angular velocity (scalar part 0).
    omegaQuat = [0; omegaBody(:)];

    % Rate of change of orientation.
    dq = 0.5 * quatMultiply(q, omegaQuat);

    % Baumgarte-style drift correction toward |q| = 1. normGain sets how quickly
    % norm error is pulled out; 1 is a mild, stable choice.
    normGain = 1;
    dq = dq + normGain * (1 - (q(:).' * q(:))) * q(:);
end
