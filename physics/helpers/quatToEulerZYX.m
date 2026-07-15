function euler = quatToEulerZYX(q)
% QUATTOEULERZYX  Roll-pitch-yaw (ZYX / 3-2-1) Euler angles from a
% scalar-first unit quaternion (body->world). For plotting/inspection only;
% the dynamics integration never uses Euler angles.
%
% INPUT
%   q : 4x1 quaternion [q0;q1;q2;q3].
%
% OUTPUT
%   euler : 1x3 [roll, pitch, yaw] (rad).

    q = q(:) / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    roll  = atan2( 2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2) );
    pitch = asin( max(-1, min(1, 2*(q0*q2 - q3*q1) )) );   % clamp: guard asin domain
    yaw   = atan2( 2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2) );

    euler = [roll, pitch, yaw];
end