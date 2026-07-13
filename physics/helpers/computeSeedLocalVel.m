function seedLocalVels = computeSeedLocalVel(seedParams, currComPos, seedRotQuat, comVel, comOmega)
% COMPUTESEEDLOCALVEL  Local flow velocity at each strip geometric centre.
%
% Computes, for one timestep, the body-frame velocity of every strip geometric
% centre (the fixed mid-chord/centroid reference, NOT the centre of pressure),
% the relative wind it sees, and the in-plane (spanwise-projected) wind used by
% the 2D strip aerodynamics.
%
% INPUTS
%   seedParams  : struct from setupSeedShapeAndMass (strip geometry + mass props).
%   currComPos  : 3x1 CoM position WITHIN THE SEED (body-datum coords, m) at
%                 the current frame. The caller interpolates this from the
%                 time-varying massParams.com_t in the main dynamics loop and
%                 passes it in. It is the reference point for the omega x r
%                 lever arms of the transport theorem.
%   seedRotQuat : 4x1 unit quaternion [q0;q1;q2;q3], scalar-first (Hamilton),
%                 representing the BODY->WORLD rotation (body axes expressed in
%                 the inertial frame). This is the orientation integrated by ode45.
%   comVel      : 3x1 CoM translational velocity in the INERTIAL frame (m/s).
%   comOmega    : 3x1 body angular velocity, expressed in BODY-frame axes (rad/s).
%
% OUTPUT (struct, all fields 3xM where M = number of strips, all body frame)
%   seedLocalVels.totalVel      : velocity of each geometric centre (body, m/s).
%   seedLocalVels.windVel       : relative wind = -totalVel (body frame, m/s).
%   seedLocalVels.projectedWind : windVel with the spanwise (body z) component
%                                 removed; rows 1-2 (x=chordwise, y=normal) are
%                                 the in-plane wind for the 2D strip aero.
%
% =========================================================================
% FRAME / QUATERNION REASONING  (read before trusting the math)
% -------------------------------------------------------------------------
% Convention: q = [q0;q1;q2;q3] is scalar-first and represents BODY->WORLD,
% i.e. a vector known in body coordinates is mapped to inertial coordinates by
%       v_inertial = R(q) * v_body.
% R(q) is built below in quatToRotm. Because R is a rotation (orthonormal),
% its inverse is its transpose, so the reverse mapping is
%       v_body = R(q)' * v_inertial.
%
% Velocity of a material point P (a strip geometric centre) on a rigid body,
% relative to the inertial frame, is the rigid-body transport relation:
%       v_P = v_CoM + omega x r_{P/CoM}                                  (1)
% This is a statement about physical vectors; it holds in any single frame as
% long as every term is expressed in that frame.
%
% We want (1) in the BODY frame. The two given pieces are in different frames:
% comVel is inertial, comOmega is already body. Rotating (1) into the body
% frame and using the identity  R'(a x b) = (R'a) x (R'b)  (valid for any
% rotation R) gives:
%       v_P^body = R' v_CoM^inertial  +  omega^body x r_{P/CoM}^body      (2)
% So we only need to rotate comVel; comOmega and the geometry are used as-is.
% This matches stripVelocityBody in the dynamics skeleton.
%
% Note on currComPos: it is the CoM location in body-datum coords at this
% frame, supplied by the caller (interpolated from massParams.com_t). It is
% the reference point subtracted from each strip geometric centre to form the
% lever arm r_{P/CoM}^body in eqn (2). It carries no inertial-position info --
% velocity is independent of absolute position -- so it is used only for the
% relative geometry, not added anywhere as a translation.
% =========================================================================

% -------------------------------------------------------------------------
% 1. ORIENTATION: build body<->world rotation from the quaternion
% -------------------------------------------------------------------------
% Normalise defensively: ode45 integration lets the quaternion norm drift, and
% a non-unit quaternion would scale (not just rotate) the velocities.
q = seedRotQuat(:) / norm(seedRotQuat);

R_body2world = quatToRotm(q);     % v_inertial = R_body2world * v_body
R_world2body = R_body2world';     % v_body     = R_world2body * v_inertial

% -------------------------------------------------------------------------
% 2. CoM VELOCITY: inertial -> body frame  (the first term of eqn (2))
% -------------------------------------------------------------------------
comVel_body = R_world2body * comVel(:);   % 3x1, body frame

% -------------------------------------------------------------------------
% 3. STRIP GEOMETRY: geometric-centre positions relative to the CoM, body frame
% -------------------------------------------------------------------------
xgc = seedParams.strips.xgc_body;   % 1xM chordwise geometric centre, body x (m)
zgc = seedParams.strips.zgc_body;   % 1xM spanwise  geometric centre, body z (m)
numStrips = numel(xgc);

% CoM location in BODY-DATUM coords, supplied by the caller for this frame.
comPos_body = currComPos(:);   % 3x1, body-datum coords (m)

% -------------------------------------------------------------------------
% 4. PER-STRIP VELOCITIES
% -------------------------------------------------------------------------
totalVel      = zeros(3, numStrips);   % geometric-centre velocity, body frame
windVel       = zeros(3, numStrips);   % relative wind = -velocity
projectedWind = zeros(3, numStrips);   % wind with spanwise (z) component removed

for i = 1 : numStrips

    % Geometric-centre position relative to the CoM, body frame (the lever arm
    % r_{P/CoM} in eqn (2)). Body y = 0 because the plate is flat (body x-z plane).
    r_gc_body = [xgc(i); 0; zgc(i)] - comPos_body;

    % Total velocity of the geometric centre in the body frame:
    %   v_P^body = comVel_body + omega^body x r_{P/CoM}^body   (eqn (2))
    v_gc_body = comVel_body + cross(comOmega(:), r_gc_body);
    totalVel(:, i) = v_gc_body;

    % Relative wind seen by the strip in still fluid = -(strip velocity).
    wind = -v_gc_body;
    windVel(:, i) = wind;

    % Projected wind for 2D strip theory: discard the spanwise (body z)
    % component, leaving the in-plane (chordwise x, normal y) flow.
    projectedWind(:, i) = [wind(1); wind(2); 0];

end

% -------------------------------------------------------------------------
% 5. PACK OUTPUT
% -------------------------------------------------------------------------
seedLocalVels.totalVel      = totalVel;        % 3xM, body frame (m/s)
seedLocalVels.windVel       = windVel;         % 3xM, body frame (m/s)
seedLocalVels.projectedWind = projectedWind;   % 3xM, body frame, z = 0 (m/s)

end   % computeSeedLocalVel





% function seedLocalVels = computeSeedLocalVel(seedParams, currComPos, seedRotQuat, comVel, comOmega)
%     %this function computes the local velocity at each seed strip center of pressure, for this
%     %timestep
%     %INPUTS: 
%     % - the seed params with all the seed shape info
%     % - the current COM position of the seed at this time
%     % - the rotation quaternion that defines the seed body frame position
%     % W.R.T the inertial frame (this comes out of ode45 integration)
%     % - the velocity of the COM (in inertial frame)
%     % - the omega of the COM(defined in the body frame axes
% 
% 
%     %Loop through each seed strip in seedParams
%         %Convert the inertial-frame COM velocity into the body frame
%         %add the contribution of the rotational velocity to each strip
%         %store this total velocity
%         %take the negative, this is the total wind direction, store this 
%         %take the total wind direction, and remove the spanwise component
%         %(Z component). This is the projected wind velocity, store it
%     %end
% end