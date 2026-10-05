function fig = plotPaperComparison(xaGrid, ybGrid, newMode, refMode, track, ttl)
% PLOTPAPERCOMPARISON  Our paper-window mode map vs Hou et al. (2025) Fig. 2a, at a glance.
%
% Three panels over their axes (x_c/a spanwise, y_c/b chordwise CoM offset):
%   1. OURS, six-part labels -- seedModePalette colours, one-letter code per cell.
%   2. OURS MAPPED TO THEIR FAMILIES -- each label translated through
%      seedModePalette().paper into AR / ST / CH / FA. Non-physical cells are dark
%      grey and labels with no paper equivalent light grey; both count as
%      disagreements. Every cell that disagrees with panel 3 has a RED outline.
%   3. THEIRS -- the digitised reference (paperModeReference), same family colours.
% Panels 2 and 3 use the family palette of reportPaperModeGrid, so the two maps
% are directly comparable by colour.
%
% INPUTS
%   xaGrid (1xNx), ybGrid (1xNy) : the grid axes
%   newMode : NyxNx string, classifySeedMode labels
%   refMode : NyxNx string, 'AR'|'ST'|'CH'|'FA'
%   track   : summarizeModeGrid output for the same grid (agreement, non-physical)
%   ttl     : extra title text (e.g. label and git hash)
% OUTPUT
%   fig : figure handle (invisible); the caller exports and closes it.

    pal = seedModePalette();
    famCol = containers.Map({'AR','ST','CH','FA'}, ...
        {[0.95 0.75 0.15], [0.85 0.25 0.25], [0.25 0.45 0.85], [0.75 0.55 0.85]});
    colNonPhys = [0.30 0.30 0.30];   colNoEquiv = [0.82 0.82 0.82];
    [nY, nX] = size(newMode);

    % ours -> family (or non-physical / no equivalent)
    fam = strings(nY, nX);  isNP = false(nY, nX);
    for iy = 1:nY
        for ix = 1:nX
            j = find(strcmp(pal.names, newMode(iy,ix)), 1);
            isNP(iy,ix) = ~pal.physical(j);
            fam(iy,ix)  = string(pal.paper{j});
        end
    end
    match = ~isNP & fam == refMode;

    fig = figure('Visible','off','Color','w','Position',[60 60 1650 520]);
    tl  = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf(['Hou et al. window -- agreement %.0f%% (%d/%d), non-physical %.0f%%' ...
                       '      %s'], 100*track.paperAgreeFrac, track.paperAgreeCount, ...
                       track.nValid, 100*track.nonPhysicalFrac, ttl), ...
          'Interpreter','none', 'FontWeight','bold');
    names = {'Ours: six-part classifier', 'Ours, mapped to their families', ...
             'Hou et al. Fig. 2a (digitised)'};
    for panel = 1:3
        ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
        for iy = 1:nY
            for ix = 1:nX
                edge = [0.35 0.35 0.35];  lw = 0.5;
                switch panel
                    case 1
                        j = find(strcmp(pal.names, newMode(iy,ix)), 1);
                        c = pal.colors(j,:);
                    case 2
                        if isNP(iy,ix);            c = colNonPhys;
                        elseif fam(iy,ix) == "";   c = colNoEquiv;
                        else;                      c = famCol(char(fam(iy,ix)));
                        end
                        if ~match(iy,ix);  edge = [0.9 0 0];  lw = 2.2;  end
                    case 3
                        c = famCol(char(refMode(iy,ix)));
                end
                scatter(ax, xaGrid(ix), ybGrid(iy), 170, 's', 'filled', ...
                        'MarkerFaceColor', c, 'MarkerEdgeColor', edge, 'LineWidth', lw);
                if panel == 1
                    text(ax, xaGrid(ix), ybGrid(iy), pal.codes(j), 'HorizontalAlignment', ...
                         'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', ...
                         [1 1 1] * (mean(c) < 0.55));
                end
            end
        end
        xlabel(ax, 'x_c/a   (spanwise CoM offset)');
        ylabel(ax, 'y_c/b   (chordwise CoM offset)');
        title(ax, names{panel});
        xlim(ax, [0 0.43]);  ylim(ax, [0 0.215]);

        % legends via dummy handles
        if panel == 1
            present = unique(newMode(:)).';
            h = gobjects(0);  lab = {};
            for nm = present
                j = find(strcmp(pal.names, nm), 1);
                h(end+1) = scatter(ax, nan, nan, 60, 's', 'filled', 'MarkerFaceColor', ...
                                   pal.colors(j,:), 'MarkerEdgeColor', [0.35 0.35 0.35]); %#ok<AGROW>
                lab{end+1} = sprintf('%s  %s', pal.codes(j), nm); %#ok<AGROW>
            end
            legend(ax, h, lab, 'Location', 'eastoutside', 'Interpreter', 'none');
        elseif panel == 2
            fams = {'AR','ST','CH','FA'};
            h = gobjects(0);  lab = {};
            for k = 1:4
                h(end+1) = scatter(ax, nan, nan, 60, 's', 'filled', 'MarkerFaceColor', ...
                                   famCol(fams{k}), 'MarkerEdgeColor', [0.35 0.35 0.35]); %#ok<AGROW>
                lab{end+1} = fams{k}; %#ok<AGROW>
            end
            h(end+1) = scatter(ax, nan, nan, 60, 's', 'filled', 'MarkerFaceColor', colNonPhys);
            lab{end+1} = 'non-physical';
            h(end+1) = scatter(ax, nan, nan, 60, 's', 'filled', 'MarkerFaceColor', colNoEquiv);
            lab{end+1} = 'no paper equivalent';
            h(end+1) = scatter(ax, nan, nan, 60, 's', 'MarkerEdgeColor', [0.9 0 0], 'LineWidth', 2);
            lab{end+1} = 'disagrees with Hou';
            legend(ax, h, lab, 'Location', 'southoutside', 'NumColumns', 4);
        end
    end
end
