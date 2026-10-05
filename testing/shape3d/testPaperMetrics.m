%% Regression: the paper-comparison metrics and classifier, on SYNTHETIC motion
% computePaperMetrics and classifyPaperMode decide every conclusion in the
% phase-5 comparison, so they are tested against trajectories whose answer is
% known by construction rather than against simulation output (which would only
% tell us the two agree with each other).
%
% HOW THE SYNTHETIC MOTION IS BUILT
% Each case prescribes a rigid-body motion analytically and then derives the state
% the metrics consume, so there is no integration and no fitting anywhere:
%
%   R(t) = Ry_world(phi(t)) * Rx_world(beta) * Rz_body(psi(t))
%     phi   revolution about the world vertical (their phi)
%     beta  fixed tilt, putting the LONG axis beta degrees below horizontal
%     psi   self-rotation about the body long axis (their psi -- the tumbling)
%
%   span axis in world = Ry(phi)*Rx(beta)*[0;0;1], whose vertical component is
%   -sin(beta) and is INDEPENDENT of psi -- i.e. tumbling does not change the cone,
%   which is the property the whole ST/AR distinction rests on.
%
%   Differentiating R gives the body angular velocity exactly:
%       omega_body = phi_dot * (R' * Yhat)  +  psi_dot * zhat
%   and since (R'*Yhat)(3) = -sin(beta),
%       omega_z = psi_dot - phi_dot*sin(beta).
%   So a PURE REVOLUTION (psi_dot = 0) at a tilt still has a constant nonzero
%   omega_z. That is the trap case A exists to catch.
%
%   Position: a helix of radius R_h about the vertical, descending at V, generated
%   as Ry(phi)*[R_h;0;0] -- by the same rotation as the attitude, so the orbit
%   turns in the same sense as the body, as it does for a real revolving seed.
%
% CASES AND WHAT EACH PROVES
%   A  pure revolution, 11.4 deg tilt (their AR case)      -> AR
%      and, critically, NOT CST: omega_z is 7.9 rad/s here purely from the
%      revolution projecting onto the tilted span axis. This is the defect that
%      writing the test found -- the classifier read omega_z directly and called
%      a clean autorotation "continuous spiral tumbling".
%   B  revolution + steady tumbling, 7 tumbles/rev,
%      38.2 deg tilt (their continuous ST case)            -> CST
%   C  the same with psi_dot reversing sign periodically
%      (their segmented ST turning points)                 -> SST
%   D  straight tilted fall, no rotation at all            -> FA
%   E  unsettled, wildly swinging attitude                 -> CH
%   F  metric accuracy: theta, phi', psi', tumbles/rev, helix radius and the
%      reversal count all recovered to the prescribed values
%   G  the sign convention matches theirs: nut-side tip below horizontal is
%      NEGATIVE theta
%   N1-N9  the six-part classifier (classifySeedMode): autorotation vs
%      spiralGlide at the 1-span radius, tumble with a level and a TILTED span
%      (no revolution either way), tight vs wide spiral tumbling, flutter,
%      end-on fall, parachute

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

L = 0.060;                                   % plate length, for R/L
mopts = struct('windowStartFrac', 0.5, 'convergeTol', 0.20, 'refLength', L);
results = {};

%% A. Pure revolution at a tilt -> AR (and must NOT read as tumbling)
[tA, xA] = synthSeedMotion(struct('beta',11.4, 'phiDot',40, 'psiDot',0, ...
                              'Rh',0.02, 'V',1.3, 'T',6));
pmA = computePaperMetrics(tA, xA, mopts);
modeA = classifyPaperMode(pmA);
wzA   = mean(abs(xA(:,13)));                 % raw |omega_z| -- the trap
okA = strcmp(modeA,'AR') && pmA.selfRotationRate < 1e-6 && wzA > 5;
results(end+1,:) = {'A pure revolution at 11.4 deg -> AR (not CST)', okA, ...
    sprintf('mode %s; raw |w_z| = %.2f rad/s but psi'' = %.1e', modeA, wzA, pmA.selfRotationRate)};

%% B. Revolution + steady tumbling -> CST
beta = 38.2;  phiDot = 6.0;  tumbTarget = 7;
psiDotB = tumbTarget * phiDot;        % psi' IS the prescribed self-rotation rate;
                                      % omega_z is what differs from it, by
                                      % phi_dot*sin(beta)
[tB, xB] = synthSeedMotion(struct('beta',beta, 'phiDot',phiDot, 'psiDot',psiDotB, ...
                              'Rh',1.25*L, 'V',1.6, 'T',12));
pmB = computePaperMetrics(tB, xB, mopts);
modeB = classifyPaperMode(pmB);
okB = strcmp(modeB,'CST') && pmB.tumbleReversals == 0;
results(end+1,:) = {'B revolution + steady tumbling -> CST', okB, ...
    sprintf('mode %s, %.2f tumbles/rev, %d reversals', modeB, pmB.tumblesPerRev, pmB.tumbleReversals)};

%% C. The same, with the tumbling direction reversing -> SST
[tC, xC] = synthSeedMotion(struct('beta',beta, 'phiDot',phiDot, 'psiDot',psiDotB, ...
                              'Rh',1.89*L, 'V',1.8, 'T',12, 'reversePeriod',3.0));
pmC = computePaperMetrics(tC, xC, mopts);
modeC = classifyPaperMode(pmC);
okC = strcmp(modeC,'SST') && pmC.tumbleReversals >= 2;
results(end+1,:) = {'C tumbling with direction reversals -> SST', okC, ...
    sprintf('mode %s, %d reversals', modeC, pmC.tumbleReversals)};

%% D. Straight tilted fall, no rotation -> FA
[tD, xD] = synthSeedMotion(struct('beta',30, 'phiDot',0, 'psiDot',0, ...
                              'Rh',0, 'V',3.0, 'T',6));
pmD = computePaperMetrics(tD, xD, mopts);
modeD = classifyPaperMode(pmD);
okD = strcmp(modeD,'FA');
results(end+1,:) = {'D tilted fall, no rotation -> FA', okD, ...
    sprintf('mode %s, theta %+.1f, V_d %.2f', modeD, pmD.spanAxisTiltDeg, pmD.descentSpeed)};

%% E. Unsettled, wildly swinging attitude -> CH
[tE, xE] = synthSeedMotion(struct('chaotic', true, 'T', 12));
pmE = computePaperMetrics(tE, xE, mopts);
modeE = classifyPaperMode(pmE);
okE = strcmp(modeE,'CH');
results(end+1,:) = {'E unsettled, swinging attitude -> CH', okE, ...
    sprintf('mode %s, tilt spread %.1f deg, converged %d', modeE, pmE.spanAxisTiltStd, pmE.converged)};

%% H. Fast rocking about the long axis that never completes a turn -> FL
% Amplitude 40 deg at 7 Hz: peak |psi'| is ~31 rad/s, well above the tumbling rate
% threshold, but the accumulated one-way rotation is only 0.11 turns. This is the
% motion our model actually produces over the left half of the paper window, and
% it must NOT be reported as FA (which means a large steady tilt).
[tH, xH] = synthSeedMotion(struct('beta',0, 'phiDot',0, 'psiDot',0, 'Rh',0, 'V',1.4, ...
                              'T',8, 'rockAmpDeg',40, 'rockHz',7));
pmH = computePaperMetrics(tH, xH, mopts);
modeH = classifyPaperMode(pmH);
okH = strcmp(modeH,'FL') && pmH.selfRotationRate > 5 && pmH.maxTurnsOneWay < 1;
results(end+1,:) = {'H fast rocking, never a full turn -> FL (not FA)', okH, ...
    sprintf('mode %s, psi'' %.1f rad/s, %.2f turns one way, %d reversals', ...
            modeH, pmH.selfRotationRate, pmH.maxTurnsOneWay, pmH.tumbleReversals)};

%% F. Metric accuracy against the prescribed values (case B)
errTheta = abs(pmB.spanAxisTiltDeg - (-beta));
errPhi   = abs(pmB.revolutionRate - phiDot) / phiDot;
errPsi   = abs(pmB.selfRotationRate - psiDotB) / abs(psiDotB);
errTumb  = abs(pmB.tumblesPerRev - tumbTarget) / tumbTarget;
errRad   = abs(pmB.helixRadiusPerL - 1.25) / 1.25;
errSd    = pmB.spanAxisTiltStd;
okF = errTheta < 1e-9 && errPhi < 1e-6 && errPsi < 1e-6 && errTumb < 1e-6 ...
      && errRad < 1e-6 && errSd < 1e-9;
results(end+1,:) = {'F metrics recover the prescribed motion', okF, ...
    sprintf('theta %.1e deg, phi'' %.1e, psi'' %.1e, tumb/rev %.1e, R/L %.1e', ...
            errTheta, errPhi, errPsi, errTumb, errRad)};

%% G. Sign convention matches theirs (nut-side tip down = negative theta)
[tG, xG] = synthSeedMotion(struct('beta',-25, 'phiDot',5, 'psiDot',0, 'Rh',0.05, 'V',1.5, 'T',6));
pmG = computePaperMetrics(tG, xG, mopts);
okG = pmB.spanAxisTiltDeg < 0 && pmG.spanAxisTiltDeg > 0 ...
      && abs(pmG.spanAxisTiltDeg - 25) < 1e-9;
results(end+1,:) = {'G theta sign: nut-side tip down is negative', okG, ...
    sprintf('beta=+38.2 -> %+.1f deg,  beta=-25 -> %+.1f deg', ...
            pmB.spanAxisTiltDeg, pmG.spanAxisTiltDeg)};

%% N. The six-part classifier (classifySeedMode) on the same synthetic motion
% Each case is built so its name is known by construction. Span length L, so the
% tight/wide boundary (1 span) sits at Rh = L. N4 is the trap that motivated
% measuring revolution from the span HEADING: a tumble about a tilted span axis
% projects onto the body's spin about the vertical, but is not a revolution.
so = struct('windowStartFrac', 0.5, 'convergeTol', 0.20, 'spanLength', L);
nCases = { ...
  'N1 revolution, radius 0.33 span -> autorotation', 'autorotation', ...
      struct('beta',40, 'phiDot',20, 'psiDot',0,  'Rh',0.33*L, 'V',1.5, 'T',8); ...
  'N2 revolution, radius 2.5 span -> spiralGlide',   'spiralGlide', ...
      struct('beta',40, 'phiDot',10, 'psiDot',0,  'Rh',2.5*L,  'V',1.5, 'T',8); ...
  'N3 tumble, level span, no revolution -> tumble',  'tumble', ...
      struct('beta',0,  'phiDot',0,  'psiDot',40, 'Rh',0,      'V',1.5, 'T',8); ...
  'N4 tumble, span tilted 25 deg, no revolution -> tumble', 'tumble', ...
      struct('beta',25, 'phiDot',0,  'psiDot',40, 'Rh',0,      'V',1.5, 'T',8); ...
  'N5 tumble + revolution within a span -> tightSpiralTumble', 'tightSpiralTumble', ...
      struct('beta',25, 'phiDot',8,  'psiDot',60, 'Rh',0.2*L,  'V',1.8, 'T',8); ...
  'N6 tumble + revolution beyond a span -> spiralTumble', 'spiralTumble', ...
      struct('beta',25, 'phiDot',8,  'psiDot',60, 'Rh',1.5*L,  'V',1.8, 'T',8); ...
  'N7 rocking, never a full turn -> flutter',        'flutter', ...
      struct('beta',0,  'phiDot',0,  'psiDot',0,  'Rh',0,      'V',1.4, 'T',8, ...
             'rockAmpDeg',40, 'rockHz',7); ...
  'N8 span 80 deg below horizontal, no rotation -> endOnFall', 'endOnFall', ...
      struct('beta',80, 'phiDot',0,  'psiDot',0,  'Rh',0,      'V',6.0, 'T',6); ...
  'N9 broadside, straight down -> parachute',        'parachute', ...
      struct('beta',0,  'phiDot',0,  'psiDot',0,  'Rh',0,      'V',1.3, 'T',6)};
for c = 1:size(nCases,1)
    [tN, xN] = synthSeedMotion(nCases{c,3});
    smN = computeSeedModeMetrics(tN, xN, so);
    [mN, pN] = classifySeedMode(smN);
    okN = strcmp(mN, nCases{c,2});
    results(end+1,:) = {nCases{c,1}, okN, sprintf('got %s [%s/%s/%s], heading rate %.1f', ...
        mN, pN.flip, pN.revolution, pN.attitude, smN.spanHeadingRate)}; %#ok<SAGROW>
end

%% Report
fprintf('\n');
allOk = true;
for k = 1:size(results,1)
    fprintf('  %-52s %-5s %s\n', results{k,1}, pf(results{k,2}), results{k,3});
    allOk = allOk && results{k,2};
end
if allOk
    fprintf('\nPASS: the paper metrics recover prescribed motion exactly, and the\n');
    fprintf('      classifier labels each synthetic mode correctly.\n');
else
    fprintf(2, '\nFAIL: see the checks above.\n');
end


% The synthetic motions come from testing/helpers/synthSeedMotion.m (promoted
% from a local function here so the mode-example gallery can share them).

function s = pf(ok);  if ok; s = 'PASS'; else; s = 'FAIL'; end;  end
