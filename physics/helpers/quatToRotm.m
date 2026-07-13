function R = quatToRotm(q)
% QUATTOROTM  Rotation matrix from a scalar-first unit quaternion (Hamilton).
%
% Convention matches the rest of the seed code: q = [q0; q1; q2; q3] with q0 the
% scalar part, representing the BODY->WORLD rotation. So a body-frame vector maps
% to the inertial frame by  v_inertial = R * v_body.
%
% INPUT
%   q : 4x1 quaternion [q0; q1; q2; q3]. Normalised internally (defensive).
%
% OUTPUT
%   R : 3x3 rotation matrix, body->world.

    q = q(:) / norm(q);                 % defensive normalisation
    w = q(1); x = q(2); y = q(3); z = q(4);

    R = [ 1 - 2*(y^2 + z^2),   2*(x*y - w*z),     2*(x*z + w*y);
          2*(x*y + w*z),       1 - 2*(x^2 + z^2), 2*(y*z - w*x);
          2*(x*z - w*y),       2*(y*z + w*x),     1 - 2*(x^2 + y^2) ];
end
