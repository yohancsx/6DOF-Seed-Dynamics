function res = testActivePathClean()
% TESTACTIVEPATHCLEAN  No active code calls archived code.
%
% Retired code (the old classifier, the planar model, duplicate scripts) is moved
% into an archive folder, which is NOT on the path. A leftover call to it would only
% fail when that line runs -- possibly deep inside a long suite run. This static
% check finds such calls up front.
%
% THE ARCHIVED NAMES are read from the archive itself: every .m file under
% testing/archive/ and archive/ (if present). Nothing is hardcoded, so archiving a
% file automatically extends the check.
%
% THE ACTIVE CODE is every .m file in the repo outside those folders and outside
% 'Olivia Code/'. Comments are ignored (a comment may name archived code), and so
% are these deliberate users of the archive, which put it on the path themselves:
%   testing/classifier/buildClassifierReviewSet.m, buildRelabelSheet.m
%       their sheets carry the old classifier's labels as comparison columns;
%   any file listed in opts.exempt (none by default).
% A file that adds an archive folder to the path itself (addpath(...'archive'...))
% is also exempt: it has made the dependency explicit.
%
% OUTPUT
%   res : struct array (.name, .pass, .detail), one entry per active file that
%         calls archived code, or a single PASS entry.

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  P = seedPaths(here);
    root = P.root;

    % --- archived names ----------------------------------------------------
    arch = [dir(fullfile(root, 'testing', 'archive', '**', '*.m')); ...
            dir(fullfile(root, 'archive', '**', '*.m'))];
    names = unique(erase({arch.name}, '.m'));

    % --- active files --------------------------------------------------------
    mf  = dir(fullfile(root, '**', '*.m'));
    rel = cellfun(@(d, n) erase(fullfile(d, n), [root filesep]), {mf.folder}, {mf.name}, ...
                  'UniformOutput', false);
    rel = strrep(rel, '\', '/');
    skip = startsWith(rel, {'testing/archive/', 'archive/', 'Olivia Code/'}) | ...
           ismember(rel, {'testing/classifier/buildClassifierReviewSet.m', ...
                          'testing/classifier/buildRelabelSheet.m', ...
                          'testing/shape3d/testActivePathClean.m'});
    rel = rel(~skip);

    res = struct('name', {}, 'pass', {}, 'detail', {});
    if isempty(names)
        res(1) = struct('name', 'active path clean', 'pass', true, 'detail', 'archive is empty');
        return
    end
    pat = ['(?<![\w.])(' strjoin(regexptranslate('escape', names), '|') ')(?=\s*\(|\s*;|\s*$|\s*,)|@(' ...
           strjoin(regexptranslate('escape', names), '|') ')\>'];
    for k = 1:numel(rel)
        txt = fileread(fullfile(root, rel{k}));
        if ~isempty(regexp(txt, 'addpath\([^\n]*archive', 'once')); continue; end
        code = stripComments(txt);
        hits = regexp(code, pat, 'match');
        if ~isempty(hits)
            hits = unique(erase(hits, '@'));
            res(end+1) = struct('name', rel{k}, 'pass', false, ...
                                'detail', ['calls archived ' strjoin(hits, ', ')]); %#ok<AGROW>
        end
    end
    if isempty(res)
        res(1) = struct('name', 'active path clean', 'pass', true, 'detail', ...
            sprintf('%d active files, %d archived names, no calls', numel(rel), numel(names)));
    end

    for k = 1:numel(res)
        if res(k).pass; tag = 'PASS'; else; tag = 'FAIL'; end
        fprintf('  %-48s %s  %s\n', res(k).name, tag, res(k).detail);
    end
    if all([res.pass])
        fprintf('PASS: no active code calls archived code.\n');
    else
        fprintf(2, 'FAIL: %d active file(s) call archived code.\n', numel(res));
    end
end


function out = stripComments(txt)
% Drop MATLAB comments: block comments, then everything from a % that is not
% inside a single- or double-quoted string. (A transpose quote after a name,
% bracket or another quote is not the start of a string.)
    txt   = regexprep(txt, '(?m)^\s*%\{.*?^\s*%\}', '');
    lines = splitlines(txt);
    for i = 1:numel(lines)
        L = lines{i};  inS = '';  cut = numel(L);
        for j = 1:numel(L)
            ch = L(j);
            if isempty(inS)
                if ch == '%'; cut = j - 1; break; end
                if ch == '"'; inS = '"';
                elseif ch == '''' && (j == 1 || isempty(regexp(L(j-1), '[\w\)\]\}''.]', 'once')))
                    inS = '''';
                end
            elseif ch == inS
                inS = '';
            end
        end
        lines{i} = L(1:cut);
    end
    out = strjoin(lines, newline);
end
