function L = crossTabLines(rowLab, colLab, ttl)
% CROSSTABLINES  A labelled count table of two label arrays, as text lines.
%
% For summaries: how often each row label co-occurs with each column label over
% the same cells (e.g. our six-part label vs the paper-vocabulary label vs the
% digitised reference map), so where two scorings disagree is visible at a glance.
%
% INPUTS
%   rowLab, colLab : same-size string/cellstr arrays of labels (one per cell).
%   ttl            : (optional) a title line.
% OUTPUT
%   L : column string array of lines (title, header, one line per row label).

    a = string(rowLab(:));   b = string(colLab(:));
    ua = unique(a);   ub = unique(b);
    w  = max([strlength(ua); 8]) + 2;
    cw = max([strlength(ub); 4]) + 2;
    L = strings(0, 1);
    if nargin >= 3 && strlength(string(ttl)) > 0; L(end+1) = "  " + string(ttl); end
    L(end+1) = "  " + pad("", w) + join(pad(ub.', cw, 'left'), "");
    for i = 1:numel(ua)
        n = arrayfun(@(j) nnz(a == ua(i) & b == ub(j)), 1:numel(ub));
        L(end+1) = "  " + pad(ua(i), w) + join(pad(string(n), cw, 'left'), ""); %#ok<AGROW>
    end
    L = L(:);
end
