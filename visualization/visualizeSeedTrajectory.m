function visualizeSeedTrajectory(t, posWorld, quatOut, opts)
% VISUALIZESEEDTRAJECTORY  Two-panel flight summary: 3D CoM trajectory (left,
% larger) and Euler angles vs time (right).
%
% AXIS CONVENTION (read before trusting the panels):
% World frame is Y-up (gravity = [0;-g;0]), but MATLAB's plot3 renders its
% THIRD argument as the vertical screen axis. visualizeSeedShape.m and
% animateSeed.m both handle this with toPlot = @(w) [w(1,:); w(3,:); w(2,:)]
% -- i.e. world Y (vertical) is routed into plot3's Z-slot, world Z
% (spanwise) into its Y-slot -- so "up" actually renders as up. The left
% panel here uses that identical remap.
%
% The SAME remap is applied to the Euler angles on the right. quatToEulerZYX
% ties roll/pitch/yaw to the x/y/z quaternion components by construction
% (roll<-q1/x, pitch<-q2/y, yaw<-q3/z), so swapping pitch and yaw mirrors the
% position swap exactly: the world-Y-associated angle (pitch) is drawn in the
% same "slot" as world Y is in the trajectory panel (3rd), and the world-Z-
% associated angle (yaw) in the same slot as world Z (2nd). This keeps the
% two panels reading consistently -- rotation about the vertical axis is
% "pitch" here, not "yaw" as it would be in a Z-down aerospace convention.
%
% INPUTS
%   t        : 1xN or Nx1 time vector (s).
%   posWorld : 3xN CoM position, inertial frame (m), rows [worldX;worldY;worldZ].
%   quatOut  : 4xN orientation quaternion history [q0;q1;q2;q3], body->world,
%              scalar-first.
%   opts     : (optional) struct of display toggles. Defaults:
%     .fig1Name   - figure name                          (default: 'Seed flight history')
%     .fig2Name   - figure name                          (default: 'Seed angle history')
%     .lineWidth - line width for both panels             (default: 1.5)
%     .view      - [az el] camera angle for the 3D panel  (default: [35 20])
%
% OUTPUT
%   none -- creates a figure with the two panels.

% -------------------------------------------------------------------------
% 0. DEFAULT OPTIONS
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = setDefault(opts, 'fig1Name',   'Seed trajectory history');
opts = setDefault(opts, 'fig2Name',   'Seed angle history');
opts = setDefault(opts, 'lineWidth', 1.5);
opts = setDefault(opts, 'view',      [35 20]);

t    = t(:);          % Nx1
numT = numel(t);

% -------------------------------------------------------------------------
% 1. EULER ANGLES from the quaternion history (raw order: roll, pitch, yaw)
% -------------------------------------------------------------------------
eulerRaw = zeros(numT, 3);
for k = 1 : numT
    eulerRaw(k, :) = quatToEulerZYX(quatOut(:, k));
end

% -------------------------------------------------------------------------
% 2. AXIS REMAP (world Y <-> world Z), applied identically to position and
%    to the pitch/yaw Euler angles -- see file header for the reasoning.
% -------------------------------------------------------------------------
posPlot   = [posWorld(1,:); posWorld(3,:); posWorld(2,:)];    % 3xN
eulerPlot = [eulerRaw(:,1),  eulerRaw(:,3),  eulerRaw(:,2)];  % Nx3: [roll, pitch, yaw]

% -------------------------------------------------------------------------
% 3. First Figure: 3D CoM trajectory
% -------------------------------------------------------------------------
figure('Name', opts.fig1Name);

plot3(posPlot(1,:), posPlot(2,:), posPlot(3,:), 'LineWidth', opts.lineWidth);
grid on; 
daspect([1 1 1])
xlabel('World X (m)'); ylabel('World Z (m)'); zlabel('World Y (m) -- up');
title('CoM trajectory (inertial frame)');
view(opts.view);

% -------------------------------------------------------------------------
% 4. Second Figure: Euler angles vs time
% -------------------------------------------------------------------------
figure('Name', opts.fig2Name);
plot(t, eulerPlot(:,1), 'LineWidth', opts.lineWidth); hold on;
plot(t, eulerPlot(:,2), 'LineWidth', opts.lineWidth);
plot(t, eulerPlot(:,3), 'LineWidth', opts.lineWidth);
grid on;
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('Roll (about X)', 'Pitch (about Z, spanwise)', 'Yaw (about Y, vertical)', ...
       'Location', 'best');
title('Euler angles');

end   % visualizeSeedTrajectory


% =========================================================================
% LOCAL: set default struct field if absent
% =========================================================================
function s = setDefault(s, field, value)
    if ~isfield(s, field)
        s.(field) = value;
    end
end
