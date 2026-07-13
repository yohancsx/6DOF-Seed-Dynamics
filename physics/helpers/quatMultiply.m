function r = quatMultiply(p, q)
% QUATMULTIPLY  Hamilton product r = p (x) q of two scalar-first quaternions.
%
% Both quaternions use the [w; x; y; z] (scalar-first) convention. The product
% composes rotations: if p and q are body->world rotations, p(x)q applies q
% first, then p.
%
% INPUTS
%   p : 4x1 quaternion [pw; px; py; pz].
%   q : 4x1 quaternion [qw; qx; qy; qz].
%
% OUTPUT
%   r : 4x1 quaternion product, scalar-first.

    pw = p(1); px = p(2); py = p(3); pz = p(4);
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);

    r = [ pw*qw - px*qx - py*qy - pz*qz;    % scalar
          pw*qx + px*qw + py*qz - pz*qy;    % x
          pw*qy - px*qz + py*qw + pz*qx;    % y
          pw*qz + px*qy - py*qx + pz*qw ];  % z
end
