function fig = plotSeedCaseTrajectory(result, visible)
% PLOTSEEDCASETRAJECTORY  One 3D CoM-trajectory figure for a single case.
%
% Uses the same world->plot remap as the rest of the codebase
% (toPlot = [X; Z; Y]) so world Y (gravity/vertical) renders as the vertical
% screen axis. Marks the start (green o) and end (red x) so the descent
% direction is unambiguous.
%
% INPUTS
%   result  : a result struct from runOneSeedCase (uses .x, .group, .label).
%   visible : 'on' or 'off' for the figure (default 'off' for batch runs).
%
% OUTPUT
%   fig : handle to the created figure.

    if nargin < 2 || isempty(visible); visible = 'off'; end

    pos     = result.x(:, 1:3).';                    % 3xN world CoM position
    posPlot = [pos(1,:); pos(3,:); pos(2,:)];        % remap: X, Z, Y-up

    fig = figure('Name', result.label, 'Color', 'w', 'Visible', visible);
    plot3(posPlot(1,:), posPlot(2,:), posPlot(3,:), 'LineWidth', 1.5);
    hold on;
    plot3(posPlot(1,1),   posPlot(2,1),   posPlot(3,1),   'go', ...
          'MarkerFaceColor', 'g', 'MarkerSize', 7);    % start
    plot3(posPlot(1,end), posPlot(2,end), posPlot(3,end), 'rx', ...
          'LineWidth', 1.5, 'MarkerSize', 9);          % end

    grid on; axis equal;
    xlabel('World X (m)'); ylabel('World Z (m)'); zlabel('World Y (m) -- up');
    title(sprintf('%s / %s', result.group, result.label), 'Interpreter', 'none');
    view(35, 20);
end
