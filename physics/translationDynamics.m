function a_inertial = translationDynamics(F_total_inertial, R, mp)
% TRANSLATIONDYNAMICS  CoM linear acceleration in the inertial frame.
%
% Solves  (M*I3 + A_inertial) * a = F_total,  where the body-frame translational
% added-mass block is rotated into the inertial frame:  A_inertial = R*A_trans*R'.
% Working in the inertial frame means no rotating-frame (Coriolis) terms appear;
% the 2D code carried those because it wrote the balance in the body frame.
%
% The added-mass rate term (-dA_inertial/dt * v) is neglected: it scales with the
% fluid density and is tiny for a seed in air.
%
% INPUTS
%   F_total_inertial : 3x1 net force in the inertial frame (aero + gravity) (N).
%   R                : 3x3 body->world rotation matrix.
%   mp               : mass-properties struct (uses mp.M, mp.A_trans).
%
% OUTPUT
%   a_inertial : 3x1 CoM linear acceleration, inertial frame (m/s^2).

    A_inertial = R * mp.A_trans * R.';          % added mass in inertial frame
    M_eff      = mp.M * eye(3) + A_inertial;     % effective mass matrix

    a_inertial = M_eff \ F_total_inertial;
end
