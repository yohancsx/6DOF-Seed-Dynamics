%% Relabelling sheet -- one dropdown column per part, pre-filled from a first pass
% Turns a labelled review set (from buildClassifierReviewSet) into a second-pass
% labelling sheet where each of the six parts of the classifier is judged in its
% own DROPDOWN column, so the same motion always gets the same name.
%
% It does NOT re-simulate: it reads the saved trajectories and the first-pass
% spreadsheet from the review folder, and writes a NEW workbook next to them
% (classifier_review_v2.xlsx, or a timestamped name if that exists). The
% first-pass sheet is never modified.
%
% THE 'review' SHEET
%   A-C    run_id, video, summary
%   D-J    YOUR PARTS, each a dropdown:
%            flip_eye        none | flutter | tumble | segmentedTumble | irregular
%            revolution_eye  none | tight | wide
%            attitude_eye    flat | steep | endOn            (span axis from horizontal)
%            path_eye        steep | glide                   (only when not revolving)
%            facing_eye      broadside | edgeOn              (only when not rotating)
%            regularity_eye  steady | wobbling | chaotic
%            physical_eye    yes | no
%   K      label_eye         dropdown of the 17 names + unsure (typing a NEW name is
%                            allowed, with a warning -- that is how a new mode starts)
%   L      name_from_parts   FORMULA: the name the classifier's precedence rules give
%                            for YOUR parts (MODE_DEFINITIONS.md section 3)
%   M      check             FORMULA: "MISMATCH" when label_eye and name_from_parts
%                            disagree -- either the label or a part needs another look
%   N-O    confidence_eye (1-3 dropdown), notes_eye
%   P-Q    your first-pass label and notes, for reference
%   R-...  THE CODE'S ANSWERS, grouped and COLLAPSED so they stay out of sight while
%          labelling (click the [+] above the columns to expand): the six-part
%          classifier's name and parts, the two older classifiers (paper metrics
%          recomputed with the 2026-10-03 sign fix), model, nut position, status
%
% PRE-FILL: each first-pass label is mapped onto ONLY the parts it actually
% states. "dive" says no rotation, steep path, edge-on -- but not whether the span
% was level, tilted or vertical, so attitude is left blank and so is the name
% (chordDive / obliqueDive / endOnFall depend on it). Nothing from the code is
% copied into your columns. Blank cells are the ones still to judge.
%
% Needs Excel (driven through actxserver) for the dropdowns, formulas and
% grouping. Section-organised, no local functions.

%% 1. Configuration  -- EDIT HERE
root      = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
reviewDir = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Classifier Review\2026-09-27_165346_3727986";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'testing','helpers'));

%% 2. Read the first pass and re-classify every saved trajectory
first = readtable(fullfile(reviewDir, 'classifier_review.xlsx'), 'Sheet', 'review', ...
                  'TextType', 'string');
txt = @(s) strtrim(lower(fillmissing(string(s), 'constant', "")));
first.label_eye = txt(first.label_eye);   first.notes_eye = txt(first.notes_eye);
first = sortrows(first, 'run_id');        % blind random order, as in the first pass
n = height(first);

ths = defaultSeedModeThresholds();  thf = defaultModeThresholds();  thp = defaultPaperThresholds();
newName = strings(n,1);  P = strings(n,6);  curName = strings(n,1);  papName = strings(n,1);
mNames = {'flipPPTurns','flipMaxTurnsOneWay','flipReversals','flipReversalCV', ...
          'spanHeadingRate','absSpanTiltDeg','spanAxisTiltStd','radiusPerSpan', ...
          'glideRatio','coneAngleDeg','spanAlongVel','descentSpeed','verticalSpinMag', ...
          'selfRotationRate','revolutionRate','converged'};
M = nan(n, numel(mNames));
for i = 1:n
    S = load(fullfile(reviewDir, 'trajectories', first.run_id(i) + ".mat"));
    if ~strcmp(S.runInfo.status, 'OK')
        newName(i) = "failed";  curName(i) = "failed";  papName(i) = "failed";  continue
    end
    so = S.cfg.metricOpts;  so.spanLength = S.cfg.spanLength;
    so.flipHysteresisTurns = ths.flipHysteresisTurns;
    sm = computeSeedModeMetrics(S.t, S.x, so);
    [newName(i), pt] = classifySeedMode(sm, ths);
    yn = ["no", "yes"];
    P(i,:) = [pt.flip, pt.revolution, pt.attitude, pt.path, pt.regularity, ...
              yn(pt.physical + 1)];
    curName(i) = classifyFlightMode(sm, thf);
    papName(i) = classifyPaperMode(sm, thp);      % with the sign fix
    for j = 1:numel(mNames);  M(i,j) = double(sm.(mNames{j}));  end
end
M(~isfinite(M)) = NaN;
fprintf('re-classified %d runs\n', n);

%% 3. Map each first-pass label onto ONLY the parts it states
% columns: flip, revolution, attitude, path, facing, regularity, physical, name
E = strings(n, 8);
for i = 1:n
    lab = first.label_eye(i);   nt = first.notes_eye(i);
    switch lab
        case "tumbling";                 E(i,:) = ["tumble","none","","","","","","tumble"];
        case "tumbling dive";            E(i,:) = ["tumble","none","","steep","","","","tumble"];
        case "tight spiral";             E(i,:) = ["tumble","tight","","","","","","tightSpiralTumble"];
        case "glide";                    E(i,:) = ["none","none","","glide","broadside","","","glide"];
        case "dive";                     E(i,:) = ["none","none","","steep","edgeOn","","",""];
        case "autorotation";             E(i,:) = ["none","tight","","","","","","autorotation"];
        case "flutter";                  E(i,:) = ["flutter","none","","","","","","flutter"];
        case "fluttering dive";          E(i,:) = ["flutter","none","","steep","","","","flutter"];
        case "fluttering spiral";        E(i,:) = ["flutter","","","","","","","flutterSpiral"];
        case {"spiral glide", "gliding spiral"}
                                         E(i,:) = ["none","wide","","","","","","spiralGlide"];
        case "spiral dive";              E(i,:) = ["none","","","","","","",""];
        otherwise                        % no first-pass label: read the notes
            if contains(nt, "chao")
                E(i,:) = ["","","","","","chaotic","no","chaotic"];
            elseif contains(nt, "weird") || contains(nt, "not usually seen")
                E(i,:) = ["","","","","","","no",""];
            elseif contains(nt, "not sure") || contains(nt, "hard to say")
                E(i,:) = ["","","","","","","","unsure"];
            end
    end
    % notes that add a part to a labelled run
    if lab ~= "" && contains(nt, "chao");          E(i,6) = "chaotic";  E(i,7) = "no";  end
    if contains(nt, "unphysical") || contains(nt, "not usually seen in real"); E(i,7) = "no"; end
end

%% 4. Write the workbook (data first; dropdowns/formulas/grouping in section 5)
outFile = fullfile(reviewDir, 'classifier_review_v2.xlsx');
if isfile(outFile)
    outFile = fullfile(reviewDir, sprintf('classifier_review_v2_%s.xlsx', ...
                       char(datetime('now','Format','yyyyMMdd_HHmmss'))));
end
blank = strings(n,1);
review = table(first.run_id, first.video, first.summary, ...
    E(:,1), E(:,2), E(:,3), E(:,4), E(:,5), E(:,6), E(:,7), E(:,8), blank, blank, ...
    nan(n,1), blank, first.label_eye, first.notes_eye, ...
    newName, P(:,1), P(:,2), P(:,3), P(:,4), P(:,5), P(:,6), curName, papName, ...
    first.model, first.nut_chordFrac, first.nut_spanFrac, first.status, ...
    'VariableNames', {'run_id','video','summary', ...
      'flip_eye','revolution_eye','attitude_eye','path_eye','facing_eye', ...
      'regularity_eye','physical_eye','label_eye','name_from_parts','check', ...
      'confidence_eye','notes_eye','label_first_pass','notes_first_pass', ...
      'label_new','flip_new','revolution_new','attitude_new','path_new', ...
      'regularity_new','physical_new','label_current','label_paper_current', ...
      'model','nut_chordFrac','nut_spanFrac','status'});
metrics = [table(first.run_id, 'VariableNames', {'run_id'}), ...
           array2table(M, 'VariableNames', mNames)];

listNames = {'flip','revolution','attitude','path','facing','regularity','physical', ...
             'label','confidence'};
lists = { ["none";"flutter";"tumble";"segmentedTumble";"irregular"], ...
          ["none";"tight";"wide"], ["flat";"steep";"endOn"], ["steep";"glide"], ...
          ["broadside";"edgeOn"], ["steady";"wobbling";"chaotic"], ["yes";"no"], ...
          ["endOnFall";"endOnSpin";"chaotic";"wobblingSpin";"edgeSlide";"chordDive"; ...
           "obliqueDive";"glide";"parachute";"autorotation";"spiralGlide";"flutter"; ...
           "flutterGlide";"flutterSpiral";"tumble";"tightSpiralTumble";"spiralTumble";"unsure"], ...
          ["1";"2";"3"] };
nMax = max(cellfun(@numel, lists));
L = strings(nMax, numel(lists));
for c = 1:numel(lists);  L(1:numel(lists{c}), c) = lists{c};  end
listTable = array2table(L, 'VariableNames', listNames);

partsHelp = table( ...
    ["flip_eye";"revolution_eye";"attitude_eye";"path_eye";"facing_eye";"regularity_eye"; ...
     "physical_eye";"label_eye";"name_from_parts / check"], ...
    ["rotation about the span: none | flutter = rocks, never a full turn | tumble = whole turns one way | segmentedTumble = whole turns, reversing periodically | irregular = whole turns, reversing irregularly"; ...
     "about the vertical: none | tight = axis within 1 span of the seed | wide = radius of 1 span or more"; ...
     "span axis from horizontal: flat < 30 deg | steep 30-70 | endOn > 70 (non-physical)"; ...
     "only if NOT revolving: steep = glide ratio <= 0.25 (nearly straight down) | glide = travels sideways"; ...
     "only if NOT rotating at all: broadside = plate face down | edgeOn = plate edge leading (a dive)"; ...
     "steady | wobbling = span axis wobbles widely | chaotic = no sustained pattern (non-physical)"; ...
     "no = a motion real seeds do not show"; ...
     "the mode name (MODE_DEFINITIONS.md section 3); type a new one if none fits"; ...
     "computed from YOUR parts by the classifier's rules; MISMATCH = label and parts disagree"], ...
    'VariableNames', {'column','meaning'});

writetable(review,    outFile, 'Sheet', 'review');
writetable(metrics,   outFile, 'Sheet', 'metrics');
writetable(partsHelp, outFile, 'Sheet', 'how_to');
writetable(listTable, outFile, 'Sheet', 'lists');
fprintf('wrote data to %s\n', outFile);

%% 5. Dropdowns, formulas, grouping (Excel via actxserver)
colL = @(k) [repmat(char(64 + floor((k-1)/26)), 1, k > 26), char(65 + mod(k-1, 26))];
vn   = review.Properties.VariableNames;
col  = @(name) colL(find(strcmp(vn, name), 1));
last = n + 1;                              % data rows 2..n+1

xl = actxserver('Excel.Application');
xl.Visible = false;   xl.DisplayAlerts = false;
try
    wb = xl.Workbooks.Open(char(outFile));
    ws = wb.Worksheets.Item('review');
    wl = wb.Worksheets.Item('lists');

    % --- dropdowns: parts STOP on invalid entries; the name only WARNS ------
    eyeCols = {'flip_eye','revolution_eye','attitude_eye','path_eye','facing_eye', ...
               'regularity_eye','physical_eye','label_eye','confidence_eye'};
    for c = 1:numel(eyeCols)
        lc  = colL(c);   nL = numel(lists{c});
        src = sprintf('=lists!$%s$2:$%s$%d', lc, lc, nL + 1);
        rg  = ws.Range(sprintf('%s2:%s%d', col(eyeCols{c}), col(eyeCols{c}), last));
        rg.Validation.Delete();
        alert = 1;                                    % xlValidAlertStop
        if strcmp(eyeCols{c}, 'label_eye'); alert = 2; end   % xlValidAlertWarning
        rg.Validation.Add(3, alert, 1, src);          % xlValidateList, ..., xlBetween
        rg.Validation.IgnoreBlank = true;
        rg.Validation.InCellDropdown = true;
    end

    % --- name_from_parts: the classifier's precedence rules on YOUR parts ---
    r = 2;
    F  = sprintf('%s%d', col('flip_eye'), r);       R = sprintf('%s%d', col('revolution_eye'), r);
    A  = sprintf('%s%d', col('attitude_eye'), r);   Pp = sprintf('%s%d', col('path_eye'), r);
    Fa = sprintf('%s%d', col('facing_eye'), r);     G = sprintf('%s%d', col('regularity_eye'), r);
    f = ['=IF(OR(' F '="",' R '="",' A '=""),"",' ...
         'IF(AND(' A '="endOn",' F '="none",' R '="none"),"endOnFall",' ...
         'IF(' A '="endOn","endOnSpin",' ...
         'IF(OR(' G '="chaotic",' F '="irregular"),"chaotic",' ...
         'IF(AND(' G '="wobbling",' R '<>"none"),"wobblingSpin",' ...
         'IF(AND(' F '="none",' R '="none"),' ...
             'IF(' Fa '="edgeOn",IF(' A '="flat","chordDive","obliqueDive"),' ...
                'IF(' Pp '="glide","glide","parachute")),' ...
         'IF(' F '="none",IF(' R '="tight","autorotation","spiralGlide"),' ...
         'IF(' F '="flutter",IF(' R '="none",IF(' Pp '="glide","flutterGlide","flutter"),"flutterSpiral"),' ...
         'IF(' R '="none","tumble",IF(' R '="tight","tightSpiralTumble","spiralTumble"))))))))))'];
    ws.Range(sprintf('%s2:%s%d', col('name_from_parts'), col('name_from_parts'), last)).Formula = f;
    K  = sprintf('%s%d', col('label_eye'), r);   Nf = sprintf('%s%d', col('name_from_parts'), r);
    ws.Range(sprintf('%s2:%s%d', col('check'), col('check'), last)).Formula = ...
        ['=IF(OR(' K '="",' Nf '=""),"",IF(' K '=' Nf ',"","MISMATCH"))'];

    % --- look: header bold, eye columns shaded, code columns grouped + collapsed
    ws.Rows.Item(1).Font.Bold = true;
    ws.Range(sprintf('%s1:%s%d', col('flip_eye'), col('notes_eye'), last)).Interior.Color = ...
        hex2dec('DDF3FF');                                   % BGR: light yellow
    ws.Range(sprintf('%s1:%s%d', col('name_from_parts'), col('check'), last)).Interior.Color = ...
        hex2dec('E6E6E6');                                   % formulas: grey
    ws.Columns.AutoFit();
    ws.Range([col('notes_eye') ':' col('notes_eye')]).ColumnWidth = 40;
    ws.Range([col('notes_first_pass') ':' col('notes_first_pass')]).ColumnWidth = 40;
    codeRange = ws.Range(sprintf('%s:%s', col('label_new'), colL(numel(vn))));
    codeRange.Columns.Group();
    ws.Outline.ShowLevels(0, 1);                            % collapse the code columns

    ws.Activate;
    xl.ActiveWindow.SplitRow = 1;   xl.ActiveWindow.SplitColumn = 1;
    xl.ActiveWindow.FreezePanes = true;
    wb.Save();   wb.Close(false);
catch ME
    try wb.Close(false); catch; end
    xl.Quit;  delete(xl);
    rethrow(ME);
end
xl.Quit;  delete(xl);

nPre  = sum(any(E ~= "", 2));   nName = sum(E(:,8) ~= "");
fprintf(['done: %s\n  %d of %d rows pre-filled from the first pass (%d with a name);\n' ...
         '  blank cells are the ones still to judge. Code columns are collapsed (click [+]).\n'], ...
        outFile, nPre, n, nName);
