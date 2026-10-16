function appendTrackingRow(file, row)
% APPENDTRACKINGROW  Add one run's row to a tracking-history CSV, merging columns.
%
% writetable(..., 'WriteMode','append') writes the new row positionally and
% cannot add columns, so a run that tracks a NEW number (a new stage, a new
% metric) would silently mis-align the history. This reads the history, adds any
% column either side is missing (NaN for numbers, "" for text), and rewrites the
% file: old rows keep their values, new columns start at this row.
%
% INPUTS
%   file : the history CSV (created if it does not exist).
%   row  : a one-row table of this run's numbers.

    if ~isfile(file)
        writetable(row, file);
        return
    end
    % Pin the header line and every type: with only a few rows, auto-detection can
    % miss the header, and it reads the timestamps as datetimes.
    hdr   = readlines(file);
    names = cellstr(strsplit(strtrim(hdr(1)), ','));
    io = delimitedTextImportOptions('NumVariables', numel(names), 'VariableNames', names, ...
                                    'Delimiter', ',', 'DataLines', [2 Inf]);
    textVars = intersect(io.VariableNames, {'timestamp', 'git', 'label', 'profile', ...
                                            'levApplicationPoint'});
    io = setvartype(io, setdiff(io.VariableNames, textVars), 'double');
    if ~isempty(textVars); io = setvartype(io, textVars, 'string'); end
    old = readtable(file, io);

    for v = setdiff(row.Properties.VariableNames, old.Properties.VariableNames)
        old.(v{1}) = blankLike(row.(v{1}), height(old));
    end
    for v = setdiff(old.Properties.VariableNames, row.Properties.VariableNames)
        row.(v{1}) = blankLike(old.(v{1}), height(row));
    end
    row = row(:, old.Properties.VariableNames);
    for v = old.Properties.VariableNames          % logical vs numeric read-back
        if islogical(row.(v{1})) && isnumeric(old.(v{1})); row.(v{1}) = double(row.(v{1})); end
        if isstring(old.(v{1})) && ~isstring(row.(v{1})); row.(v{1}) = string(row.(v{1})); end
    end
    writetable([old; row], file);
end

function b = blankLike(x, n)
    if isstring(x) || iscellstr(x) || ischar(x)
        b = repmat("", n, 1);
    else
        b = nan(n, 1);
    end
end
