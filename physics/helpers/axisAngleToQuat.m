function q = axisAngleToQuat(axis, angleRad)
% AXISANGLETOQUAT  Scalar-first unit quaternion for a rotation about an axis.
%
% Returns q = [cos(a/2); sin(a/2)*axisUnit], the body->world convention used by
% the rest of the seed code (quatToRotm / quatToEulerZYX). Convenient for
% building single-axis initial orientations (e.g. an initial pitch or roll tilt)
% without composing full Euler angles.
%
% INPUTS
%   axis     : 3x1 rotation axis (need not be unit; normalised internally). A
%              (near) zero-length axis returns the identity quaternion.
%   angleRad : rotation angle about that axis (rad).
%
% OUTPUT
%   q : 4x1 unit quaternion [q0; q1; q2; q3], scalar-first.

    axis = axis(:);
    n = norm(axis);
    if n < eps
        q = [1; 0; 0; 0];              % undefined axis -> identity
        return
    end
    u = axis / n;                       % unit rotation axis
    q = [cos(angleRad/2); sin(angleRad/2) * u];
    q = q / norm(q);                    % defensive re-normalisation
end
