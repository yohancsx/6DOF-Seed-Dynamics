%% Explore the LEV (leading-edge vortex) term -- plotted from the IMPLEMENTED code
% Evaluates the shape3d model's LEV vortex-lift term on its own, straight from
% physics3d/computeLEVForce and physics3d/levPlanformConstants, so the figures show
% exactly what the model does when enableLEV is on. Five figures:
%   1. the vortex-lift increment vs angle of attack, against the existing APW law
%   2. the Rossby gate vs Rossby number, with where THIS seed's strips actually sit
%   3. the gated increment over the whole (alpha, Ro) plane
%   4. the two application points ('colocated' vs 'forward')
%   5. the gate vs REVOLUTION RATE -- what the default (kinematic) gate reads
% Figures are written OUTSIDE the repo (generated binaries are kept out of git).
%
% Every equation, its source, and every modelling choice is documented at its
% point of use in computeLEVForce and levPlanformConstants. In brief:
%
%   C_L,v = K_v sin^2(a) cos(a)          Rezgui, Arroyo & Theunissen (2020),
%                                        Aeronautical J. 124(1278):1236-1261, eq. (3)
%   K_v   = K_p - K_p^2 K_i               ibid. eq. (4); K_p Helmbold (1942),
%                                        K_i = 1/(pi AR) Prandtl elliptic
%   G     = 1/(1 + (Ro/Ro_crit)^p)        OUR construction (see below)
%   Ro    = |v_ip|/(Omega c)  'kinematic' (DEFAULT) Rossby's own ratio;
%                                        Omega = |omega x s_hat|, the strip's
%                                        revolution rate about the CoM
%         = r/c               'geometric' the pure-revolution special case, kept
%                                        as a toggle (lev.rossbyDefinition)
%   x_app = l_cp*c            'colocated' Snyder & Lamar (1972), NASA TN D-6994
%   x_app = lambda_v (c/2) cos(a) 'forward'  no direct citation
%
% ON THE ROSSBY NUMBER -- verified against Lentink & Dickinson (2009), J. Exp. Biol.
% 212:2705-2719. They define Ro = Rg/c (radius of gyration over MEAN chord) for
% their robot experiments, observing a stable LEV on a revolving wing at Ro = 2.9,
% an unstable one on a translating wing (Ro = inf), and force coefficients that
% change with Ro over 2.9 -> 3.6 -> 4.4. For the survey of real wings they switch
% to the TIP radius, Ro = R/c (their Eq. 6, equal to the single-wing aspect ratio),
% and find wings cluster near 3 -- which is Rg/c ~ 1.5 for insects. They give NO
% critical Rossby number. So Ro_crit = 3 is an anchor at the right scale, not a
% published threshold, and p has no literature value. Rezgui et al. apply the vortex
% lift everywhere, ungated (equivalent to G = 1).
%
% WHY THE DEFAULT IS THE KINEMATIC FORM. Both of their definitions are RADIUS over
% chord, which is the special case of Rossby's ratio U/(Omega L) for a wing in PURE
% REVOLUTION, where U = Omega*R. What they are actually scaling is the Coriolis and
% centripetal accelerations against the fluid's convective acceleration, and that
% general ratio is |v|/(Omega c). Read geometrically, r/c says a strip 3 chords out
% has a stable vortex whether or not the seed is revolving -- so a parachuting,
% gliding or fluttering seed gets full vortex lift, which it cannot hold. Read
% kinematically, the same strip gets Ro -> inf, G -> 0 when the seed stops
% revolving, and the two agree EXACTLY when it does revolve (testNewPhysicsTerms
% L3b/L3c). Figures 2-4 below are the same for either reading -- G(Ro) is one
% curve; only the meaning of Ro changes. Figure 5 shows what that means in Omega.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

outDir = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Outputs\LEV Exploration";
if ~exist(outDir,'dir'); mkdir(outDir); end

% --- The test seed, built the normal way so the constants are the model's own --
S = 0.050;  c = 0.015;  th = 0.002;  rho = 1.225;
b.seedShape = polyshape([-S/2,S/2,S/2,-S/2],[-c/2,-c/2,c/2,c/2]);
b.seedDensity = 65*th;  b.seedThickness = th;  b.numStrips = 10;
b.tSamples = 0;  b.nutPos_t = [0;0;0];  b.nutMass_t = 75e-6;
sp  = buildSeedParams(b, struct('rhoFluid',rho,'g',9.81,'shapeModel','shape3d'));
lev = sp.lev;
dz  = sp.strips.dz(1);

aDeg = linspace(0, 90, 361);  a = deg2rad(aDeg);

%% K_v for this seed, and its sensitivity to what "the wing" means
fprintf('\n=== LEV constants the model uses for this seed (seedParams.lev) ===\n');
fprintf('  AR = %.3f (total span / mean chord)  K_p = %.3f  K_i = %.4f  K_v = %.3f\n', ...
        lev.AR, lev.Kp, lev.Ki, lev.Kv);
fprintf('  Ro_crit = %g   p = %g   lambda_v = %g   rossbyDefinition = %s\n', ...
        lev.RoCrit, lev.p, lev.lambdaV, lev.rossbyDefinition);
fprintf('\n  Other readings of the aspect ratio (override with cfg.lev.AR):\n');
fprintf('  %-38s %-7s %-9s %-9s %-9s\n','aspect-ratio choice','AR','K_p','K_i','K_v');
ARs = {'one half-span blade from a centred CoM', (S/2)/c ; ...
       'whole seed (the model default)',         S/c     ; ...
       'Rezgui et al.''s real samara',            4.38    };
Kv = zeros(size(ARs,1),1);
for k = 1:size(ARs,1)
    L = levPlanformConstants(ARs{k,2});  Kv(k) = L.Kv;
    fprintf('  %-38s %-7.2f %-9.3f %-9.4f %-9.3f\n', ARs{k,1}, L.AR, L.Kp, L.Ki, L.Kv);
end
Kp2D = 5.2;  Ki = 1/(pi*(S/2)/c);
fprintf(['  INCONSISTENT mix (2D APW slope CL1 = 5.2 with a 3D K_i, AR = %.2f):' ...
         ' K_v = %.3f  <- why K_p and K_i must share one wing model\n'], (S/2)/c, Kp2D - Kp2D^2*Ki);

% Evaluate the IMPLEMENTED term: unit in-plane speed, flow at angle a.
% Figures 1-4 sweep Ro directly, so they drive the term through the GEOMETRIC
% definition (Ro = r/c at unit speed); G(Ro) is the same function either way, and
% this keeps "ungated" (r = 0 -> G = 1) available as a reference curve. Figure 5
% uses the model's actual default, the kinematic definition.
levGeo  = lev;  levGeo.rossbyDefinition = 'geometric';
levWith = @(L, KvVal) setfield(L, 'Kv', KvVal); %#ok<SFLD>
CTlev = @(L, aa, r) arrayfun(@(q) getfield(nthOut(3, @computeLEVForce, ...
            cos(q), sin(q), q, c, dz, r, 0, 0, L, 'colocated', rho), 'CT_lev'), aa);

co   = computeAeroCoeffs(a, []);
CT0  = co.CT;   lcp0 = co.l_cp_frac;

%% Figure 1 -- the increment vs alpha (ungated: r = 0 -> G = 1)
f1 = figure('Name','LEV vs alpha','Color','w','Position',[60 60 820 520]);
hold on; grid on;
patch([0 25 25 0],[-0.2 -0.2 3.2 3.2],[0.9 0.95 0.9],'EdgeColor','none','FaceAlpha',0.6, ...
      'DisplayName','Rezgui et al. validated range (0-25 deg)');
plot(aDeg, CT0, 'k-', 'LineWidth',2, 'DisplayName','existing APW  C_T');
cols = lines(3);
for k = 1:numel(Kv)
    d = CTlev(levWith(levGeo, Kv(k)), a, 0);
    plot(aDeg, d, '--', 'Color',cols(k,:), 'LineWidth',1.4, ...
         'DisplayName',sprintf('vortex increment, K_v=%.2f (AR %.2f)',Kv(k),ARs{k,2}));
    plot(aDeg, CT0 + d, '-', 'Color',cols(k,:), 'LineWidth',1.4, ...
         'DisplayName',sprintf('APW + vortex, K_v=%.2f',Kv(k)));
end
xline(14,':','stall 14\circ','HandleVisibility','off','LabelVerticalAlignment','middle');
xline(rad2deg(atan(sqrt(2))),':','vortex peak 54.7\circ','HandleVisibility','off', ...
      'LabelVerticalAlignment','middle');
xlabel('angle of attack \alpha (deg)'); ylabel('lift coefficient');
title('LEV vortex-lift increment, from computeLEVForce (ungated, G = 1)');
legend('Location','northwest','FontSize',8); ylim([-0.2 3.2]);
drawnow; exportgraphics(f1, fullfile(outDir,'LEV_1_vs_alpha.png'),'Resolution',140);

%% Figure 2 -- the Rossby gate, and where this seed's strips sit
Ro = linspace(0, 8, 400);  ps = [2 4 8];
f2 = figure('Name','Rossby gate','Color','w','Position',[80 80 820 460]);
hold on; grid on;
zS = sp.strips.zgc_body;
RoCentred = abs(zS)/c;          RoOffset = abs(zS - 0.06)/c;
patch([min(RoCentred) max(RoCentred) max(RoCentred) min(RoCentred)],[0 0 1.05 1.05], ...
      [0.85 0.93 1],'EdgeColor','none','FaceAlpha',0.7,'DisplayName','strips, centred CoM (collapse cases)');
patch([min(RoOffset) max(RoOffset) max(RoOffset) min(RoOffset)],[0 0 1.05 1.05], ...
      [1 0.9 0.85],'EdgeColor','none','FaceAlpha',0.7,'DisplayName','strips, nut at 1.2S (autorotation)');
for p = ps
    Lp = levGeo;  Lp.p = p;
    G  = arrayfun(@(r) getfield(nthOut(3, @computeLEVForce, 1, 1, pi/4, c, dz, r*c, 0, 0, Lp, ...
                  'colocated', rho), 'G'), Ro);
    plot(Ro, G, 'LineWidth',1.8, 'DisplayName',sprintf('p = %d%s', p, repmat(' (default)',1,p==lev.p)));
end
xline(lev.RoCrit,'k--','Ro_{crit} = 3 (anchor, not a published threshold)','HandleVisibility','off', ...
      'LabelVerticalAlignment','middle');
xlabel('local Rossby number  Ro_i   (kinematic |v_i|/(\Omega_i c_i), or geometric r_i/c_i)');
ylabel('gate  G(Ro)');
title({'Rossby gate: how much of the vortex lift a strip is allowed', ...
       'one curve for either definition -- the shaded bands are the GEOMETRIC reading'});
legend('Location','northeast','FontSize',8); ylim([0 1.05]);
drawnow; exportgraphics(f2, fullfile(outDir,'LEV_2_rossby_gate.png'),'Resolution',140);

%% Figure 3 -- gated increment over (alpha, Ro), at the model defaults
dCT = zeros(numel(Ro), numel(a));
for iR = 1:numel(Ro)
    dCT(iR,:) = CTlev(levGeo, a, Ro(iR)*c);
end
f3 = figure('Name','Gated increment','Color','w','Position',[100 100 820 520]);
imagesc(aDeg, Ro, dCT); set(gca,'YDir','normal'); colorbar;
hold on;
yline(max(RoCentred),'w--','centred-CoM strips below','LabelHorizontalAlignment','left');
yline(min(RoOffset),'w:','offset-CoM strips above','LabelHorizontalAlignment','left');
xline(25,'w-','validated \leq 25\circ');
xlabel('angle of attack \alpha (deg)'); ylabel('local Rossby number Ro');
title(sprintf('Gated vortex-lift increment G(Ro) K_v sin^2\\alpha cos\\alpha   (K_v = %.2f, p = %d)', ...
      lev.Kv, lev.p));
drawnow; exportgraphics(f3, fullfile(outDir,'LEV_3_alpha_Ro_map.png'),'Resolution',140);

%% Figure 4 -- the two application points
% Lift-weighted chordwise centre of pressure from mid-chord, in chords (+ toward
% the leading edge). Illustrative: it weights the lift components only, which is
% what the vortex term changes.
Lv   = CTlev(levGeo, a, 0);
xFwd = arrayfun(@(q) nthOut(2, @computeLEVForce, cos(q), sin(q), q, c, dz, 0, 0, 0, levGeo, ...
                'forward', rho), a) / c;
cpA  = (CT0.*lcp0 + Lv.*lcp0) ./ (CT0 + Lv);      % 'colocated': unchanged by construction
cpB  = (CT0.*lcp0 + Lv.*xFwd) ./ (CT0 + Lv);      % 'forward'
f4 = figure('Name','Application point','Color','w','Position',[120 120 820 500]);
hold on; grid on;
patch([0 25 25 0],[-0.05 -0.05 0.55 0.55],[0.9 0.95 0.9],'EdgeColor','none','FaceAlpha',0.6, ...
      'DisplayName','validated lift range (0-25 deg)');
plot(aDeg, lcp0, 'k-', 'LineWidth',2, 'DisplayName','existing APW l_{cp}(\alpha)');
plot(aDeg, cpA, '--', 'LineWidth',1.6, 'DisplayName','''colocated'' (default; Snyder & Lamar 1972)');
plot(aDeg, cpB, '-', 'LineWidth',1.6, 'DisplayName',sprintf('''forward'', \\lambda_v = %.1f (no citation)', lev.lambdaV));
plot(aDeg, xFwd, ':', 'LineWidth',1.2, 'DisplayName','''forward'' vortex application point');
yline(0.5,'k:','leading edge','HandleVisibility','off','LabelHorizontalAlignment','left');
xlabel('angle of attack \alpha (deg)'); ylabel('CoP from mid-chord  (fraction of chord, + toward LE)');
title({'levApplicationPoint: where the vortex force acts', ...
       '''colocated'' adds lift only; ''forward'' also shifts the centre of pressure'});
legend('Location','northeast','FontSize',8); ylim([-0.05 0.55]);
drawnow; exportgraphics(f4, fullfile(outDir,'LEV_4_application_point.png'),'Resolution',140);

%% Figure 5 -- the DEFAULT (kinematic) gate: how fast must the seed revolve?
% Same gate, plotted against the quantity it actually reads. At a given in-plane
% speed v the gate is half open at Omega* = v / (Ro_crit * c), so a seed that falls
% faster needs to revolve proportionally faster to hold its vortex. Driven through
% computeLEVForce with the model's own default definition.
Om    = linspace(0, 250, 500);
vSet  = [0.5 1 2 4 8];
f5 = figure('Name','Kinematic gate','Color','w','Position',[140 140 860 500]);
hold on; grid on;
cmap = parula(numel(vSet)+1);
for k = 1:numel(vSet)
    Gk = arrayfun(@(w) getfield(nthOut(3, @computeLEVForce, vSet(k), 0, 0, c, dz, 0, w, 0, ...
                  lev, 'colocated', rho), 'G'), Om);
    plot(Om, Gk, 'LineWidth',1.8, 'Color',cmap(k,:), ...
         'DisplayName',sprintf('|v| = %.1f m/s   (half open at \\Omega = %.0f rad/s)', ...
                               vSet(k), vSet(k)/(lev.RoCrit*c)));
end
% Revolution rates this model actually reaches, measured as median |omega x s_hat|
% over the second half of the reference trajectories (LEV off).
xline(17.8, 'k--', 'autorotation, 17.8', 'HandleVisibility','off','LabelVerticalAlignment','bottom');
xline(39.8, 'k:',  'tight spiral, 39.8', 'HandleVisibility','off','LabelVerticalAlignment','bottom');
xlabel('strip revolution rate  \Omega = |\omega \times s| (rad/s)');
ylabel('gate  G');
title({sprintf('Kinematic gate (default): G vs revolution rate, c = %.0f mm, Ro_{crit} = %g, p = %g', ...
       c*1e3, lev.RoCrit, lev.p), ...
       'dashed/dotted: rates the model reaches in the reference cases (measured)'});
legend('Location','southeast','FontSize',8); ylim([0 1.05]);
drawnow; exportgraphics(f5, fullfile(outDir,'LEV_5_kinematic_gate.png'),'Resolution',140);

%% Numbers worth reading off
iA = @(d) find(aDeg >= d, 1);
fprintf('\n=== At representative angles (model defaults, ungated) ===\n');
fprintf('%-8s %-10s %-10s %-9s %-13s %-13s\n','alpha','APW C_T','vortex','ratio','CoP colocated','CoP forward');
for d = [10 20 25 35 45 55 70]
    i = iA(d);
    fprintf('%-8d %-10.3f %-10.3f %-9.2f %-+13.3f %-+13.3f\n', d, CT0(i), Lv(i), ...
            Lv(i)/max(CT0(i),eps), cpA(i), cpB(i));
end
fprintf('\nGate at the defaults (Ro_crit = %g, p = %g):\n', lev.RoCrit, lev.p);
for RoQ = [0.2 1.0 1.7 3.0 4.1 5.7]
    fprintf('  Ro = %.1f -> G = %.3f\n', RoQ, 1/(1+(RoQ/lev.RoCrit)^lev.p));
end
fprintf(['\nWhat those Ro mean under the DEFAULT (kinematic) definition, c = %.0f mm:\n' ...
         '  the gate is half open at Omega* = |v|/(Ro_crit*c), i.e.\n'], c*1e3);
for vQ = [0.5 1 2 4 8]
    fprintf('  |v| = %4.1f m/s -> Omega* = %6.1f rad/s (%.1f rev/s)\n', ...
            vQ, vQ/(lev.RoCrit*c), vQ/(lev.RoCrit*c)/(2*pi));
end
fprintf(['  Measured in this model (median |omega x s|): autorotation 17.8 rad/s,\n' ...
         '  tight spiral 39.8 rad/s, and ~0 for the parachuting/gliding/fluttering\n' ...
         '  cases -- which the geometric definition would have given full vortex lift.\n']);
fprintf('\nFigures written to %s\n', outDir);

function out = nthOut(n, fn, varargin)
% Return the n-th output of fn.
    outs = cell(1, n);
    [outs{:}] = fn(varargin{:});
    out = outs{n};
end
