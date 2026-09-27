function reportPaperModeGrid(gridFile)
% REPORTPAPERMODEGRID  Read a saved paper-window mode grid and report it.
%
% Split out of runPaperModeGrid so the expensive part (the integrations) is never
% at risk from a bug in the cheap part (printing and plotting): the grid saves to
% disk first, then this runs. Re-run it alone on a saved .mat to change the
% reporting without repeating ~25 minutes of ode45.
%
% Produces:
%   1. ASCII maps -- ours, and the digitised Fig. 2a, in the same layout
%   2. agreement with the digitised reference, overall and per published region
%   3. a cross-tab of our seven-mode vocabulary against their four
%   4. a two-panel figure in their axes and palette
%
% INPUT  gridFile : path to the .mat written by runPaperModeGrid.

    S = load(gridFile);
    xaGrid = S.xaGrid;  ybGrid = S.ybGrid;
    paperMode = S.paperMode;  ourMode = S.ourMode;  refMode = S.refMode;
    nY = numel(ybGrid);  nX = numel(xaGrid);

    ourCode = containers.Map({'AR','CST','SST','CH','FA','FL','failed'}, ...
                             {'A','c','s','X','F','L','?'});
    refCode = containers.Map({'AR','ST','CH','FA'}, {'A','T','X','F'});

    %% 1. ASCII maps
    fprintf('\n=== OUR map  (rows = y_c/b, descending; cols = x_c/a) ===\n');
    fprintf('    A=AR  c=CST  s=SST  X=CH  F=FA  L=FL(flutter, not their mode)  ?=failed\n');
    for iy = nY:-1:1
        row = blanks(nX);
        for ix = 1:nX;  row(ix) = ourCode(char(paperMode(iy,ix)));  end
        fprintf('  y=%.3f | %s\n', ybGrid(iy), row);
    end
    fprintf('\n=== THEIR map, digitised  (T = ST, both kinds) ===\n');
    for iy = nY:-1:1
        row = blanks(nX);
        for ix = 1:nX;  row(ix) = refCode(char(refMode(iy,ix)));  end
        fprintf('  y=%.3f | %s\n', ybGrid(iy), row);
    end
    fprintf('           x_c/a from %.2f to %.2f\n', xaGrid(1), xaGrid(end));

    %% 2. Agreement (their ST covers both of our ST kinds)
    ourFamily = paperMode;
    ourFamily(ourFamily=="CST" | ourFamily=="SST") = "ST";
    valid = ourFamily ~= "failed";
    agree = sum(ourFamily(valid) == refMode(valid));
    fprintf('\n=== Agreement with the digitised Fig. 2a ===\n');
    fprintf('  %d / %d cells (%.0f%%)\n', agree, nnz(valid), 100*agree/max(nnz(valid),1));
    fprintf('\n  per published region:\n');
    for cat = ["ST","CH","AR","FA"]
        sel = valid & (refMode == cat);
        if nnz(sel) == 0;  continue;  end
        got = ourFamily(sel);
        u = unique(got);
        parts = cell(1, numel(u));
        for j = 1:numel(u);  parts{j} = sprintf('%s=%d', u(j), sum(got == u(j)));  end
        fprintf('    they say %-3s (%3d cells) -> we say %s\n', cat, nnz(sel), strjoin(parts, ', '));
    end

    %% 3. Cross-tab, our vocabulary vs theirs
    fprintf('\n=== Cross-tab: our label (rows) vs paper label (cols) ===\n');
    oursU  = unique(ourMode(valid));   paperU = unique(paperMode(valid));
    fprintf('%-14s', '');
    for b = 1:numel(paperU);  fprintf('%-6s', paperU(b));  end
    fprintf('\n');
    for a = 1:numel(oursU)
        fprintf('%-14s', oursU(a));
        for b = 1:numel(paperU)
            fprintf('%-6d', sum(ourMode(valid)==oursU(a) & paperMode(valid)==paperU(b)));
        end
        fprintf('\n');
    end

    %% 4. Figure
    palette = containers.Map({'AR','CST','SST','ST','CH','FA','FL','failed'}, ...
        {[0.95 0.75 0.15], [0.85 0.25 0.25], [0.85 0.45 0.45], [0.85 0.25 0.25], ...
         [0.25 0.45 0.85], [0.75 0.55 0.85], [0.35 0.75 0.45], [0 0 0]});
    f = figure('Name','Paper mode grid','Color','w','Position',[80 80 1080 460]);
    for panel = 1:2
        subplot(1,2,panel);  hold on;  box on;
        if panel == 1
            src = paperMode;  ttl = 'This model';
        else
            src = refMode;    ttl = 'Hou et al. Fig. 2a (digitised)';
        end
        for iy = 1:nY
            for ix = 1:nX
                scatter(xaGrid(ix), ybGrid(iy), 130, 'filled', 's', ...
                        'MarkerFaceColor', palette(char(src(iy,ix))), ...
                        'MarkerEdgeColor', [0.3 0.3 0.3]);
            end
        end
        xlabel('x_c/a   (spanwise CoM offset)');
        ylabel('y_c/b   (chordwise CoM offset)');
        title(ttl);  xlim([0 0.43]);  ylim([0 0.215]);
    end
    sgtitle(sprintf(['Flight mode vs CoM position -- yellow AR, red ST, blue CH, ' ...
                     'purple FA   (R%d M%d L%d)'], S.cfg.enableAddedMassRate, ...
                     S.cfg.enableAddedMassMoment, S.cfg.enableLEV));
    drawnow;
    exportgraphics(f, fullfile(S.outDir, 'paper_mode_grid.png'), 'Resolution', 150);

    %% 5. Speeds
    fprintf('\nDescent speed over the grid: %.2f to %.2f m/s (scale %.2f -> %.2f-%.2f x)\n', ...
            min(S.Vd(:)), max(S.Vd(:)), S.info.Vscale, ...
            min(S.Vd(:))/S.info.Vscale, max(S.Vd(:))/S.info.Vscale);
    fprintf('Figure written to %s\n', S.outDir);
end
