%% Explore the candidate LEV (leading-edge vortex) term BEFORE implementing it
% A self-contained PROTOTYPE: nothing here is wired into physics3d/. It evaluates
% the candidate vortex-lift term on its own so the modelling choices can be seen
% and argued over first. Four figures:
%   1. the vortex-lift increment vs angle of attack, against the existing APW law
%   2. the Rossby gate vs Rossby number, with where THIS seed's strips actually sit
%   3. the gated increment over the whole (alpha, Ro) plane
%   4. the centre-of-pressure question -- where the vortex force acts, which the
%      literature does NOT settle (see "lambda_v" below)
%
% =========================================================================
% THE FORM, AND WHERE IT COMES FROM
% =========================================================================
% Rezgui, Arroyo & Theunissen (2020), Aeronautical J. 124(1278):1236-1261,
% doi:10.1017/aer.2020.25 -- adapting Polhamus (1966, NASA TN D-3767; 1971,
% J. Aircraft 8(4):193-199, doi:10.2514/3.44254) to a samara blade section:
%
%   (1)  C_L,p = K_p sin(a) cos^2(a)                       potential lift
%   (3)  C_L,v = K_v [cos(a)/cos(Lambda)] sin^2(a)         vortex lift
%   (4)  K_v   = K_p - K_p^2 K_i
%   (5)  K_i   = d C_Di,p / d (C_L,p^2)                    induced-drag factor
%   (7)  C_l(a) = C_L,p + C_L,v   taken directly as the 2D sectional value
%
% Lambda is sweep (0 here). K_p and K_i come from lifting-surface theory for the
% planform (they used the Tornado VLM). KEY POINT: in this formulation K_v is NOT
% a free tuning constant -- it is derived from the planform. Only the ADDITION,
% eq. (3), is new relative to this model: the existing APW law already supplies a
% potential-flow lift, so adding eq. (1) as well would double-count it.
%
% Rezgui et al. validate the lift curve against Azuma & Yasuda's samara data for
% alpha in 0-25 deg only. The vortex term peaks at atan(sqrt(2)) = 54.7 deg, so
% anything past ~25 deg is extrapolation. They do NOT gate by Rossby number and
% do NOT give a centre of pressure for the vortex force -- both are additions
% made here and are flagged as such below.
%
% =========================================================================
% DEFINITIONS FOR THIS MODEL  (each is a modelling choice -- alternatives noted)
% =========================================================================
% r        Radius of strip i from the seed's rotation axis.
%          OUR CHOICE: spanwise distance from the strip's geometric centre to the
%          seed CoM, along the local span axis:  r_i = |s_i . (p_i - c_CoM)|
%          (flat seed: |z_i - c_z|). Geometric, per-strip, time-invariant for a
%          fixed CoM, tracks a moving one.
%          ALTERNATIVES: Lentink & Dickinson (2009, JEB 212:2705-2719) use ONE
%          radius per wing -- the second-moment radius of gyration measured from
%          the rotation axis, R2 = sqrt(int r^2 c dr / int c dr). Or: the
%          perpendicular distance to the INSTANTANEOUS spin axis, d x omega_hat.
%          KNOWN WEAKNESS of our choice: it is purely geometric, so it grants a
%          "stable LEV" even when the seed is not actually revolving (e.g. pure
%          fluttering). A kinematic gate would fix that -- see the open questions.
% c        Local strip chord c_i (0.015 m for every strip of the rectangular test
%          seed). Lentink & Dickinson use the mean chord; for a tapered wing the
%          local and mean chords differ, and so would the gate.
% Ro       LOCAL Rossby number Ro_i = r_i / c_i. Note this is not the same
%          quantity as Lentink & Dickinson's GLOBAL Ro = R2/c_mean: ours varies
%          along the span, theirs is one number per wing.
% Ro_crit  3. Lentink & Dickinson: Ro ~ 3 is the convergent value across insects,
%          seeds and birds; Ro = 2.9 gives a compact stable spiral LEV, Ro = inf
%          (a purely translating wing) an unstable one.
% p        Steepness of the roll-off G(Ro) = 1 / (1 + (Ro/Ro_crit)^p).
%          NO LITERATURE VALUE. Lentink & Dickinson report a threshold, not a
%          transfer function. p = 2 is gentle (G = 0.20 at 2*Ro_crit), p = 8 is
%          nearly a step. Purely a modelling choice -- figure 2 shows several.
% lambda_v Chordwise application point of the vortex force, as a fraction of the
%          half-chord from the strip centre toward the LEADING edge:
%               x_v = lambda_v * (c/2) * cos(alpha)
%          The cos(alpha) points it upwind and sends it to zero at broadside.
%          THE LITERATURE CONFLICTS WITH THE MOTIVATION FOR THIS TERM. Snyder &
%          Lamar (1972, NASA TN D-6994) find the vortex-lift load distribution on
%          low-AR delta wings is "similar in shape to that of the potential-flow
%          longitudinal loading" -- i.e. the vortex force acts at about the SAME
%          centre of pressure, not forward of it. So two options are plotted:
%            A  co-located with the existing APW centre of pressure  (no shift)
%            B  forward at lambda_v = 0.5                            (~0.28c aft
%               of the leading edge; physically motivated for an UNSWEPT wing,
%               but with no direct citation)
% K_v      Derived, per eq. (4), from K_p and K_i. Estimated here with a
%          finite-wing (Helmbold) lift slope and elliptic induced drag,
%          K_i = 1/(pi*AR). K_p and K_i MUST come from the same model: mixing the
%          2D APW slope (CL1 = 5.2) with a 3D K_i is inconsistent, and at low AR
%          it drives K_v to ~0 spuriously (printed below). AR itself is a choice:
%          the whole seed S/c, or one blade measured from the CoM.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics','aero'));

% Figures go OUTSIDE the repo (generated binaries are kept out of git).
outDir = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\LEV Exploration";
if ~exist(outDir,'dir'); mkdir(outDir); end

S = 0.050;  c = 0.015;
aDeg = linspace(0, 90, 361);  a = deg2rad(aDeg);

%% K_v from the planform -- and its sensitivity to the modelling choices
fprintf('\n=== K_v = K_p - K_p^2 K_i  (Rezgui et al. 2020, eq. 4) ===\n');
fprintf('%-34s %-7s %-9s %-9s %-9s\n','aspect-ratio choice','AR','K_p','K_i','K_v');
ARs   = {'one blade, centred CoM  (S/2)/c', (S/2)/c ; ...
         'whole seed               S/c',    S/c     ; ...
         'real samara (Rezgui)',            4.38    };
a0 = 2*pi;
Kv = zeros(size(ARs,1),1);
for k = 1:size(ARs,1)
    AR = ARs{k,2};  x = a0/(pi*AR);
    Kp = a0 / (sqrt(1 + x^2) + x);          % Helmbold finite-wing lift slope
    Ki = 1/(pi*AR);                         % elliptic induced-drag factor
    Kv(k) = Kp - Kp^2*Ki;
    fprintf('%-34s %-7.2f %-9.3f %-9.4f %-9.3f\n', ARs{k,1}, AR, Kp, Ki, Kv(k));
end
Kp2D = 5.2;  Ki = 1/(pi*(S/2)/c);
fprintf(['  INCONSISTENT mix (2D APW slope CL1=5.2 with a 3D K_i, AR=%.2f):' ...
         ' K_v = %.3f  <- spurious, do not use\n'], (S/2)/c, Kp2D - Kp2D^2*Ki);
KvRef = Kv(2);    % reference for the maps: whole-seed AR

%% Existing APW law, for comparison
co   = computeAeroCoeffs(a, []);
CT0  = co.CT;   lcp0 = co.l_cp_frac;
vort = @(Kv, a) Kv .* sin(a).^2 .* cos(a);      % eq. (3), Lambda = 0

%% Figure 1 -- the increment vs alpha
f1 = figure('Name','LEV vs alpha','Color','w','Position',[60 60 820 520]);
hold on; grid on;
patch([0 25 25 0],[-0.2 -0.2 3.2 3.2],[0.9 0.95 0.9],'EdgeColor','none','FaceAlpha',0.6, ...
      'DisplayName','Rezgui et al. validated range (0-25 deg)');
plot(aDeg, CT0, 'k-', 'LineWidth',2, 'DisplayName','existing APW  C_T');
cols = lines(3);
for k = 1:numel(Kv)
    plot(aDeg, vort(Kv(k),a), '--', 'Color',cols(k,:), 'LineWidth',1.4, ...
         'DisplayName',sprintf('vortex increment, K_v=%.2f (AR %.2f)',Kv(k),ARs{k,2}));
    plot(aDeg, CT0 + vort(Kv(k),a), '-', 'Color',cols(k,:), 'LineWidth',1.4, ...
         'DisplayName',sprintf('APW + vortex, K_v=%.2f',Kv(k)));
end
xline(14,':','stall 14\circ','HandleVisibility','off','LabelVerticalAlignment','middle');
xline(rad2deg(atan(sqrt(2))),':','vortex peak 54.7\circ','HandleVisibility','off', ...
      'LabelVerticalAlignment','middle');
xlabel('angle of attack \alpha (deg)'); ylabel('lift coefficient');
title('Candidate LEV vortex-lift increment (G = 1, i.e. ungated)');
legend('Location','northwest','FontSize',8); ylim([-0.2 3.2]);
exportgraphics(f1, fullfile(outDir,'LEV_1_vs_alpha.png'),'Resolution',140);

%% Figure 2 -- the Rossby gate, and where this seed's strips sit
RoCrit = 3;  Ro = linspace(0, 8, 400);  ps = [2 4 8];
f2 = figure('Name','Rossby gate','Color','w','Position',[80 80 820 460]);
hold on; grid on;
% per-strip Ro ranges for the two extreme configurations (10 strips)
zS = linspace(-S/2, S/2, 10);
RoCentred = abs(zS)/c;          RoOffset = abs(zS + 0.06)/c;
patch([min(RoCentred) max(RoCentred) max(RoCentred) min(RoCentred)],[0 0 1.05 1.05], ...
      [0.85 0.93 1],'EdgeColor','none','FaceAlpha',0.7,'DisplayName','strips, centred CoM (collapse cases)');
patch([min(RoOffset) max(RoOffset) max(RoOffset) min(RoOffset)],[0 0 1.05 1.05], ...
      [1 0.9 0.85],'EdgeColor','none','FaceAlpha',0.7,'DisplayName','strips, nut at 1.2S (autorotation)');
for p = ps
    plot(Ro, 1./(1+(Ro/RoCrit).^p), 'LineWidth',1.8, 'DisplayName',sprintf('p = %d',p));
end
xline(RoCrit,'k--','Ro_{crit} = 3','HandleVisibility','off');
xlabel('local Rossby number  Ro_i = r_i / c_i'); ylabel('gate  G(Ro)');
title('Rossby gate: how much of the vortex lift a strip is allowed');
legend('Location','northeast','FontSize',8); ylim([0 1.05]);
exportgraphics(f2, fullfile(outDir,'LEV_2_rossby_gate.png'),'Resolution',140);

%% Figure 3 -- gated increment over (alpha, Ro)
pRef = 4;
[AA, RR] = meshgrid(a, Ro);
dCT = (1./(1+(RR/RoCrit).^pRef)) .* vort(KvRef, AA);
f3 = figure('Name','Gated increment','Color','w','Position',[100 100 820 520]);
imagesc(aDeg, Ro, dCT); set(gca,'YDir','normal'); colorbar;
hold on;
yline(max(RoCentred),'w--','centred-CoM strips below','LabelHorizontalAlignment','left');
yline(min(RoOffset),'w:','offset-CoM strips above','LabelHorizontalAlignment','left');
xline(25,'w-','validated \leq 25\circ');
xlabel('angle of attack \alpha (deg)'); ylabel('local Rossby number Ro');
title(sprintf('Gated vortex-lift increment  \\DeltaC_T = G(Ro)K_v sin^2\\alpha cos\\alpha   (K_v=%.2f, p=%d)', KvRef, pRef));
exportgraphics(f3, fullfile(outDir,'LEV_3_alpha_Ro_map.png'),'Resolution',140);

%% Figure 4 -- the centre-of-pressure question
% Lift-weighted chordwise CoP, measured from mid-chord in units of chord
% (+ = toward the leading edge). Illustrative: it weights the lift components
% only, which is what the vortex term changes.
lamV  = 0.5;
xAPW  = lcp0;                                   % existing APW CoP (fraction of c)
xVort = lamV * 0.5 * cos(a);                    % option B location (fraction of c)
Lv    = vort(KvRef, a);
cpA   = (CT0.*xAPW + Lv.*xAPW ) ./ (CT0 + Lv);  % option A: co-located -> unchanged
cpB   = (CT0.*xAPW + Lv.*xVort) ./ (CT0 + Lv);  % option B: forward
f4 = figure('Name','Centre of pressure','Color','w','Position',[120 120 820 500]);
hold on; grid on;
patch([0 25 25 0],[-0.05 -0.05 0.55 0.55],[0.9 0.95 0.9],'EdgeColor','none','FaceAlpha',0.6, ...
      'DisplayName','validated lift range (0-25 deg)');
plot(aDeg, xAPW, 'k-', 'LineWidth',2, 'DisplayName','existing APW l_{cp}(\alpha)');
plot(aDeg, cpA, '--', 'LineWidth',1.6, 'DisplayName','option A: vortex co-located (Snyder & Lamar 1972)');
plot(aDeg, cpB, '-', 'LineWidth',1.6, 'DisplayName',sprintf('option B: vortex forward, \\lambda_v = %.1f', lamV));
plot(aDeg, xVort, ':', 'LineWidth',1.2, 'DisplayName','option B vortex application point');
yline(0.5,'k:','leading edge','HandleVisibility','off','LabelHorizontalAlignment','left');
xlabel('angle of attack \alpha (deg)'); ylabel('CoP from mid-chord  (fraction of chord, + toward LE)');
title({'Where does the vortex force act?', ...
       'This choice decides whether LEV touches the pitch balance at all'});
legend('Location','northeast','FontSize',8); ylim([-0.05 0.55]);
exportgraphics(f4, fullfile(outDir,'LEV_4_centre_of_pressure.png'),'Resolution',140);

%% Numbers worth reading off
iA = @(d) find(aDeg >= d, 1);
fprintf('\n=== At representative angles (K_v = %.2f, ungated) ===\n', KvRef);
fprintf('%-8s %-10s %-10s %-9s %-11s %-11s\n','alpha','APW C_T','vortex','ratio','CoP A','CoP B');
for d = [10 20 25 35 45 55 70]
    i = iA(d);
    fprintf('%-8d %-10.3f %-10.3f %-9.2f %-+11.3f %-+11.3f\n', d, CT0(i), Lv(i), ...
            Lv(i)/max(CT0(i),eps), cpA(i), cpB(i));
end
fprintf('\nGate G(Ro) at p = %d:\n', pRef);
for RoQ = [0.2 1.0 1.7 3.0 4.1 5.7]
    fprintf('  Ro = %.1f -> G = %.3f\n', RoQ, 1/(1+(RoQ/RoCrit)^pRef));
end
fprintf('\nFigures written to %s\n', outDir);
