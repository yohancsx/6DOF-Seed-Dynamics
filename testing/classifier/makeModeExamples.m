%% Mode examples -- one idealized animation per flight-mode label
% A visual dictionary of the 17 mode names of the six-part classifier
% (testing/classifier/MODE_DEFINITIONS.md, section 3), for anyone who has never
% seen the modes before: one short video + still per label, and an INDEX.md that
% explains each one.
%
% THE EXAMPLES ARE SYNTHETIC, NOT MODEL OUTPUT. Each is a prescribed, textbook
% version of its mode (testing/helpers/synthSeedMotion), so every label gets a
% clean example -- including ones our model never produces. Rates are chosen to be
% WATCHABLE (revolutions ~1-3 per second, tumbles ~4 turns per second), not to
% match any particular seed.
%
% SELF-CHECK: before anything is rendered, every example is run through the real
% classifier (computeSeedModeMetrics + classifySeedMode). An example that does not
% classify as its own label is reported and NOT rendered, so the gallery can never
% show a mode the code would call something else. Set cfg.render = false to run
% only the check.
%
% Videos use animateModeTrajectory in BLIND mode (the older classifier's mode
% colours would contradict the new names) with the ground-track panel; the title
% carries the mode name, its one-line definition, and NON-PHYSICAL where it applies.
%
% Output is a regenerable folder OUTSIDE the repo; reruns overwrite it.
% Section-organised, no local functions.

%% 1. Configuration  -- EDIT HERE
root   = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics";
outDir = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\Mode Examples";
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'visualization'), ...
        fullfile(root,'testing','helpers'));

cfg.render        = true;    % false -> self-check only
cfg.T             = 8;       % s of synthetic motion per example
cfg.animWindow    = 3;       % s animated, from the end
cfg.animFps       = 60;      % frames per simulated second
cfg.playbackSpeed = 0.25;    % video plays at animFps*this = 15 fps
% The seed drawn in the zoomed panel: our working seed (span sets the 1-span
% tight/wide boundary, so radii below are in units of it)
S = 0.050;   c = 0.015;

%% 2. The examples, in the precedence order of MODE_DEFINITIONS section 3
% {label, one-line definition, physical?, synthSeedMotion parameters}
% Radii are in spans of the drawn seed (S); the tight/wide boundary is 1 span.
E = {
 'endOnFall',    'falls span-first (span > 70 deg from horizontal), not turning', false, ...
     struct('beta',85,'phiDot',0,'psiDot',0,'Rh',0,'V',6.0,'T',cfg.T)
 'endOnSpin',    'falls span-first while turning', false, ...
     struct('beta',80,'phiDot',12,'psiDot',0,'Rh',0.3*S,'V',5.0,'T',cfg.T)
 'chaotic',      'no sustained pattern; the attitude swings irregularly', false, ...
     struct('chaotic',true,'T',12)
 'wobblingSpin', 'revolves, but the span axis wobbles widely', false, ...
     struct('beta',15,'phiDot',15,'psiDot',0,'Rh',0.3*S,'V',1.5,'T',cfg.T, ...
            'betaWobbleDeg',15,'betaWobbleHz',0.7)
 'edgeSlide',    'slides along its own level span, edge-on', false, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',1.0,'T',cfg.T,'psi0',90, ...
            'vHoriz',[0;0;3])
 'chordDive',    'no rotation; dives chord edge first, span level', true, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',6.0,'T',cfg.T,'psi0',90)
 'obliqueDive',  'no rotation; dives edge-on with the span tilted', true, ...
     struct('beta',45,'phiDot',0,'psiDot',0,'Rh',0,'V',6.0,'T',cfg.T,'psi0',90)
 'glide',        'no rotation; shallow attitude, travels sideways (glide ratio > 0.25)', true, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',1.0,'T',cfg.T,'psi0',-10, ...
            'vHoriz',[2;0;0])
 'parachute',    'no rotation; falls broadside, nearly straight down', true, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',1.3,'T',cfg.T)
 'autorotation', 'revolves about an axis within 1 span of the seed, no flip', true, ...
     struct('beta',15,'phiDot',18,'psiDot',0,'Rh',0.3*S,'V',1.2,'T',cfg.T)
 'spiralGlide',  'circles on a radius of 1 span or more, no flip', true, ...
     struct('beta',20,'phiDot',8,'psiDot',0,'Rh',2.5*S,'V',1.5,'T',cfg.T)
 'flutter',      'rocks about the span (never a full turn), falls nearly straight down', true, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',1.4,'T',cfg.T, ...
            'rockAmpDeg',35,'rockHz',2)
 'flutterGlide', 'rocks about the span while travelling sideways', true, ...
     struct('beta',0,'phiDot',0,'psiDot',0,'Rh',0,'V',1.4,'T',cfg.T, ...
            'rockAmpDeg',35,'rockHz',2,'vHoriz',[1.2;0;0])
 'flutterSpiral','rocks about the span while revolving', true, ...
     struct('beta',20,'phiDot',8,'psiDot',0,'Rh',2*S,'V',1.5,'T',cfg.T, ...
            'rockAmpDeg',30,'rockHz',2)
 'tumble',       'whole end-over-end turns about the span, no revolution', true, ...
     struct('beta',0,'phiDot',0,'psiDot',25,'Rh',0,'V',1.4,'T',cfg.T, ...
            'vHoriz',[0.8;0;0])
 'tightSpiralTumble', 'tumbles while revolving within 1 span', true, ...
     struct('beta',25,'phiDot',8,'psiDot',25,'Rh',0.3*S,'V',1.6,'T',cfg.T)
 'spiralTumble', 'tumbles while circling on a radius of 1 span or more', true, ...
     struct('beta',25,'phiDot',8,'psiDot',25,'Rh',2*S,'V',1.6,'T',cfg.T)
};
nE = size(E,1);

%% 3. Self-check: each example must classify as its own label
th = defaultSeedModeThresholds();
so = struct('windowStartFrac',0.5, 'convergeTol',0.20, 'spanLength',S, ...
            'flipHysteresisTurns', th.flipHysteresisTurns);
T = cell(nE,1);  X = cell(nE,1);  got = strings(nE,1);  partsStr = strings(nE,1);
ok = false(nE,1);  physGot = false(nE,1);
fprintf('\n%-3s %-18s %-18s %-40s %s\n', '#', 'intended', 'classified as', 'parts', 'ok');
for i = 1:nE
    [T{i}, X{i}] = synthSeedMotion(E{i,4});
    sm = computeSeedModeMetrics(T{i}, X{i}, so);
    [got(i), p] = classifySeedMode(sm, th);
    partsStr(i) = sprintf('%s / %s / %s / %s / %s', p.flip, p.revolution, p.attitude, ...
                          p.path, p.regularity);
    physGot(i) = p.physical;
    ok(i) = strcmp(got(i), E{i,1}) && (p.physical == E{i,3});
    fprintf('%-3d %-18s %-18s %-40s %s\n', i, E{i,1}, got(i), partsStr(i), ...
            string(ok(i)));
end
fprintf('%d of %d examples classify as their own label\n', nnz(ok), nE);
if ~all(ok)
    fprintf(2, 'NOT rendered (fix these examples or the classifier first): %s\n', ...
            strjoin(string(E(~ok,1)).', ', '));
end
if ~cfg.render; return; end

%% 4. Render the ones that passed
if ~exist(outDir, 'dir'); mkdir(outDir); end
hs = S/2;  hc = c/2;
bsp.seedShape = polyshape([-hs, hs, hs, -hs], [-hc, -hc, hc, hc]);
bsp.seedDensity = 65*0.002;  bsp.seedThickness = 0.002;  bsp.numStrips = 10;
bsp.tSamples = 0;  bsp.nutPos_t = [0;0;0];  bsp.nutMass_t = 75e-6;
sp = buildSeedParams(bsp, struct('rhoFluid',1.225,'g',9.81));   % geometry for drawing only
files = strings(nE,1);
for i = find(ok).'
    stem = sprintf('%02d_%s', i, E{i,1});
    t = T{i};  x = X{i};
    sel = t >= t(end) - cfg.animWindow;
    ttl = sprintf('%s  --  %s', E{i,1}, E{i,2});
    if ~E{i,3}; ttl = [ttl '   (NON-PHYSICAL)']; end
    animateModeTrajectory(t(sel), x(sel,:), struct( ...
        'videoFile', char(fullfile(outDir, [stem '.mp4'])), 'fps', cfg.animFps, ...
        'playbackSpeed', cfg.playbackSpeed, 'title', ttl, 'seedParams', sp, ...
        'showSeedVels', false, 'blind', true, 'topPanel', 'groundTrack'));
    exportgraphics(gcf, fullfile(outDir, [stem '.png']), 'Resolution', 110);
    close(gcf);
    files(i) = stem;
    fprintf('  rendered %s\n', stem);
end

%% 5. INDEX.md -- the reader's entry point
fid = fopen(fullfile(outDir, 'INDEX.md'), 'w');
fprintf(fid, '# Flight-mode examples\n\n');
fprintf(fid, ['One idealized example per mode name of the six-part classifier. The ' ...
              'definitions, in words and math, are in `testing/classifier/MODE_DEFINITIONS.md` ' ...
              '(section 3); the boundaries live in `testing/helpers/defaultSeedModeThresholds.m`.\n\n']);
fprintf(fid, ['**These are synthetic, not simulations.** Each is a prescribed textbook ' ...
              'version of its mode (`synthSeedMotion`), slowed to watchable rates. Every one ' ...
              'was checked to classify as its own label before it was rendered.\n\n']);
fprintf(fid, ['**Reading a video:** left, the 3D path of the centre of mass; top right, ' ...
              'the ground track seen from above (circling = revolving); bottom right, a ' ...
              'close-up of the seed at its true attitude. Videos play at %gx real time.\n\n'], ...
        cfg.playbackSpeed);
fprintf(fid, 'Generated %s by `testing/classifier/makeModeExamples.m`.\n\n', ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm')));
fprintf(fid, '| # | mode | what it looks like | physical | parts (flip / revolution / attitude / path / regularity) | files |\n');
fprintf(fid, '|---|---|---|---|---|---|\n');
for i = 1:nE
    if files(i) == ""
        fl = '*not rendered (failed self-check)*';
    else
        fl = sprintf('[video](%s.mp4) · [still](%s.png)', files(i), files(i));
    end
    phys = 'yes';  if ~E{i,3}; phys = '**no**'; end
    fprintf(fid, '| %d | `%s` | %s | %s | %s | %s |\n', i, E{i,1}, E{i,2}, phys, ...
            partsStr(i), fl);
end
fprintf(fid, ['\nNon-physical modes are motions real seeds do not show; the model can ' ...
              'produce them (see `RESEARCH_LOG.md`, phases 6-7).\n']);
fclose(fid);
fprintf('\nWritten to %s (%d videos + INDEX.md)\n', outDir, nnz(files ~= ""));
