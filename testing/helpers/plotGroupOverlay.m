function fig = plotGroupOverlay(groupResults, groupName, visible)
% PLOTGROUPOVERLAY  All trajectories of one sweep group on a single 3D axis.
%
% Overlays every successful case in a group, colour-graded in sweep order with a
% legend, so the effect of the swept parameter is visible at a glance (this is
% the figure that lets you "make generalizations"). Same world->plot remap as
% plotSeedCaseTrajectory. Non-OK cases are skipped.
%
% INPUTS
%   groupResults : struct array of result structs, all sharing one .group.
%   groupName    : the group name (title / figure name).
%   visible      : 'on' or 'off' (default 'off').
%
% OUTPUT
%   fig : handle to the created figure.

    if nargin < 3 || isempty(visible); visible = 'off'; end

    fig = figure('Name', [groupName ' overlay'], 'Color', 'w', 'Visible', visible);
    hold on;

    n    = numel(groupResults);
    cmap = turbo(max(n, 2));            % colour-grade in sweep order
    handles = gobjects(0);
    labels  = strings(0);

    for i = 1:n
        r = groupResults(i);
        if ~strcmp(r.status, 'OK') || isempty(r.x)
            continue                    % skip failed / empty
        end
        pos     = r.x(:, 1:3).';
        posPlot = [pos(1,:); pos(3,:); pos(2,:)];   % X, Z, Y-up
        h = plot3(posPlot(1,:), posPlot(2,:), posPlot(3,:), ...
                  'LineWidth', 1.5, 'Color', cmap(i,:));
        % Small end marker so overlapping trajectories are distinguishable.
        plot3(posPlot(1,end), posPlot(2,end), posPlot(3,end), ...
              'x', 'Color', cmap(i,:), 'LineWidth', 1.5, 'MarkerSize', 8);
        handles(end+1) = h;             %#ok<AGROW>
        labels(end+1)  = string(r.label); %#ok<AGROW>
    end

    grid on; 
    %axis equal;
    xlabel('World X (m)'); ylabel('World Z (m)'); zlabel('World Y (m) -- up');
    title([groupName ' -- trajectory sweep'], 'Interpreter', 'none');
    view(35, 20);
    if ~isempty(handles)
        legend(handles, labels, 'Interpreter', 'none', 'Location', 'bestoutside');
    end
end
