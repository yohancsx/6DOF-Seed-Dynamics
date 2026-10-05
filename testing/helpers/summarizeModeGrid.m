function s = summarizeModeGrid(labels, absTilt, revolving, refFamily)
% SUMMARIZEMODEGRID  The tracked numbers for one mode grid labelled by classifySeedMode.
%
% The three numbers every grid run reports, so physics changes can be compared run
% to run (roadmap: re-baseline with the fixed classifier):
%   1. NON-PHYSICAL FRACTION -- share of valid cells in a non-physical mode
%      (endOnFall, endOnSpin, chaotic, wobblingSpin, edgeSlide). Should go to 0.
%   2. SPAN TILT OF REVOLVERS -- distribution of mean |theta| (span axis from
%      horizontal) over cells that revolve (revolution part tight or wide). A real
%      samara autorotates nearly flat (Hou et al. AR: theta = -11.4 deg); steep
%      revolvers are the clearest remaining defect of the shape3d model.
%   3. PAPER AGREEMENT -- only when refFamily is given (the Hou et al. window):
%      share of valid cells whose mode maps to the digitised family
%      (seedModePalette .paper). EVERY non-physical label counts as a
%      disagreement, as MODE_DEFINITIONS.md documents; chaotic cells where the
%      paper says CH are reported separately rather than silently scored.
%
% INPUTS (same shape, any shape)
%   labels    : string array of classifySeedMode names ('failed' allowed)
%   absTilt   : mean |theta| per cell (deg), NaN where failed
%   revolving : logical, revolution part ~= 'none'
%   refFamily : (optional) string array 'AR'|'ST'|'CH'|'FA' per cell, or []
%
% OUTPUT s : struct with the numbers below and .lines (string array, ready to
%            print or append to a summary file).

    if nargin < 4; refFamily = []; end
    pal = seedModePalette();
    labels = string(labels(:));  absTilt = absTilt(:);  revolving = logical(revolving(:));

    s.nCells  = numel(labels);
    valid     = labels ~= "failed";
    s.nFailed = nnz(~valid);
    s.nValid  = nnz(valid);

    % --- census in palette order ---------------------------------------------
    counts = cellfun(@(nm) nnz(labels == nm), pal.names);
    s.census = table(string(pal.names).', counts.', 'VariableNames', {'mode','count'});

    % --- 1. non-physical fraction -------------------------------------------
    nonPhysNames = string(pal.names(~pal.physical & ~strcmp(pal.names, 'failed')));
    isNonPhys = valid & ismember(labels, nonPhysNames);
    s.nonPhysicalCount = nnz(isNonPhys);
    s.nonPhysicalFrac  = s.nonPhysicalCount / max(s.nValid, 1);

    % --- 2. span tilt of revolvers --------------------------------------------
    rv = valid & revolving & isfinite(absTilt);
    tl = absTilt(rv);
    s.revolverCount = nnz(rv);
    s.revolverFrac  = s.revolverCount / max(s.nValid, 1);
    s.tiltEdges     = 0:10:90;
    if isempty(tl)
        [s.tiltMean, s.tiltMedian, s.tiltStd, s.tiltP10, s.tiltP90, s.flatRevolverFrac] = deal(NaN);
        s.tiltHist = zeros(1, numel(s.tiltEdges)-1);
    else
        st = sort(tl);   pc = @(p) st(max(1, min(numel(st), round(p*numel(st)))));
        s.tiltMean = mean(tl);   s.tiltMedian = median(tl);   s.tiltStd = std(tl);
        s.tiltP10  = pc(0.10);   s.tiltP90 = pc(0.90);
        s.flatRevolverFrac = mean(tl < 30);       % 'flat' attitude, samara-like
        s.tiltHist = histcounts(tl, s.tiltEdges);
    end

    % --- 3. paper agreement ---------------------------------------------------
    s.hasPaper = ~isempty(refFamily);
    if s.hasPaper
        ref = string(refFamily(:));
        fam = strings(size(labels));
        for k = 1:numel(pal.names)
            fam(labels == pal.names{k}) = pal.paper{k};
        end
        match = valid & ~isNonPhys & (fam == ref);
        s.paperAgreeCount = nnz(match);
        s.paperAgreeFrac  = s.paperAgreeCount / max(s.nValid, 1);
        s.chaoticInCH     = nnz(valid & labels == "chaotic" & ref == "CH");
        cats = ["AR","ST","CH","FA"];
        s.paperByFamily = strings(numel(cats),1);
        for c = 1:numel(cats)
            sel = valid & ref == cats(c);
            if ~any(sel); continue; end
            [u,~,j] = unique(labels(sel));
            cnt = accumarray(j, 1);
            [~, o] = sort(cnt, 'descend');
            s.paperByFamily(c) = sprintf('they say %-2s (%3d): ours %s', cats(c), nnz(sel), ...
                strjoin(compose('%s=%d', u(o), cnt(o)), ', '));
        end
    end

    % --- text ------------------------------------------------------------------
    L = strings(0,1);
    L(end+1) = sprintf('  cells %d (valid %d, failed %d)', s.nCells, s.nValid, s.nFailed);
    if s.hasPaper
        L(end+1) = sprintf('  PAPER AGREEMENT      %5.1f%%  (%d/%d; non-physical counted as disagreement; chaotic in CH cells: %d)', ...
                           100*s.paperAgreeFrac, s.paperAgreeCount, s.nValid, s.chaoticInCH);
    end
    L(end+1) = sprintf('  NON-PHYSICAL         %5.1f%%  (%d/%d)', 100*s.nonPhysicalFrac, ...
                       s.nonPhysicalCount, s.nValid);
    L(end+1) = sprintf(['  REVOLVER SPAN TILT   %d revolvers (%.0f%% of valid): |theta| mean %.1f, median %.1f, ' ...
                        'p10-p90 %.1f-%.1f deg; flat (<30 deg) %.0f%%'], s.revolverCount, ...
                       100*s.revolverFrac, s.tiltMean, s.tiltMedian, s.tiltP10, s.tiltP90, ...
                       100*s.flatRevolverFrac);
    L(end+1) = "    |theta| histogram (deg): " + strjoin(compose('%d-%d:%d', ...
        s.tiltEdges(1:end-1).', s.tiltEdges(2:end).', s.tiltHist.'), '  ');
    present = s.census.count > 0;
    L(end+1) = "  census: " + strjoin(compose('%s=%d', s.census.mode(present), ...
                                       s.census.count(present)), '  ');
    if s.hasPaper
        for c = 1:numel(s.paperByFamily)
            if s.paperByFamily(c) ~= ""; L(end+1) = "    " + s.paperByFamily(c); end
        end
    end
    s.lines = L(:);    % column (growing from empty with end+1 yields a row)
end
