%% Regression: the phase-4 physics terms (added-mass rate, edge drag, LEV)
% Locks the three terms added to the shape3d model in phase 4. Defaults: edge drag
% ON, added-mass rate OFF, LEV OFF. Run after any change to rigidBody6DOF,
% computeEdgeDrag, computeLEVForce, levPlanformConstants or seed6DOFODE3D.
%
%   A. ADDED-MASS RATE (rigidBody6DOF, enableAddedMassRate) -- ON by default for
%      shape3d since Sep 2026; the frozen planar builder keeps it OFF
%      A1  the implementation equals the formula  w x (A_i v) - A_i (w x v)
%      A2  the formula equals a finite difference of R*A_b*R' along a REAL
%          ode45 trajectory -- the physics, not just the algebra
%      A3  it enters the solve with the right sign:
%          a_on - a_off = -(M*I + A_i) \ (Adot_i v)
%      A4  it vanishes identically with no rotation, and when switched off
%      A5  PROVENANCE: with it ON, the inertial acceleration equals Kirchhoff's
%          BODY-frame form  M_b vdot_b = F_b - w x (M_b v_b)  -- the form the APW
%          reference model uses (minimal_imp lines 208-209). With it OFF it does
%          not. i.e. this term is exactly what the 2D reference already contains.
%      A6  it is OFF by default in both builders
%   M. ADDED-MASS (MUNK) MOMENT (rigidBody6DOF, enableAddedMassMoment) -- also ON
%      by default for shape3d, OFF for planar
%      M1  the implementation equals  v x (A v)
%      M2  its 2D limit equals APW eq. (6.3)'s (m11 - m22) vx' vy' -- the term the
%          reference model has and this one had no counterpart for
%      M3  it enters the rotational solve with the right sign
%      M4  it vanishes along a principal axis, and when switched off
%      M5  THE PHYSICAL CLAIM: it destabilises an edge-on fall, turning the plate
%          broadside-on to its own motion
%      M6  it is OFF by default in both builders
%   B. EDGE DRAG (computeEdgeDrag, enableEdgeDrag)
%      B1-B6  zero/quadratic/odd/opposing/magnitude, and in-RHS for a pure slide
%   L. LEV VORTEX LIFT (computeLEVForce, levPlanformConstants, enableLEV)
%      L1  planform constants: the test seed's values; Helmbold -> 2*pi as AR -> inf
%          and -> pi*AR/2 as AR -> 0; Kv > 0 across aspect ratios
%      L2  the coefficient vanishes at alpha = 0 and +/-90 deg, equals
%          Kv sin^2|cos| in magnitude, and has the APW lift's sign everywhere
%      L3  the Rossby gate, 'geometric': G(0) = 1, G(RoCrit) = 1/2, monotone
%      L3b the Rossby gate, 'kinematic' (the default): it REDUCES EXACTLY to the
%          geometric gate under pure revolution (v_ip = Omega*r), and its limits are
%          well-posed -- Omega -> 0 shuts it, v_ip -> 0 opens it, both zero is
%          finite rather than 0/0
%      L3c the same two claims inside the RHS: a translating, non-rotating seed gets
%          NO vortex lift under 'kinematic' (it gets a full one under 'geometric' --
%          the defect this switch fixes), and a twisted seed in pure revolution
%          gives bit-for-bit the same gate under either definition
%      L4  the force is lift: perpendicular to the strip's in-plane velocity
%      L5  application points: 'colocated' = l_cp*c; 'forward' = lambda_v(c/2)cos(a),
%          pointing upwind and zero at broadside
%      L6  IN THE RHS: off -> bit-identical to the RHS without it; on -> the net
%          force changes by exactly sum(F_lev) and the torque by exactly
%          sum(r_app x F_lev); the two application points give the same force
%          and different torques
%      L7  defaults (off, 'colocated', 'kinematic'), the definition toggle is
%          settable through cfg.lev, and cfg.lev overrides keep Kv consistent
%   C. PLANAR IS UNTOUCHED: rate OFF, no edge-drag or LEV switch at all.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

S = 0.050;  c = 0.015;  th = 0.002;  rho = 1.225;
cfg  = struct('rhoFluid',rho,'g',9.81);
cfg3 = cfg;  cfg3.shapeModel = 'shape3d';
cfgR = cfg3; cfgR.enableAddedMassRate = true;          % rate term switched ON
b.seedShape = polyshape([-S/2,S/2,S/2,-S/2],[-c/2,-c/2,c/2,c/2]);
b.seedDensity = 65*th;  b.seedThickness = th;  b.numStrips = 10;
b.tSamples = 0;  b.nutPos_t = [c;0;1.2*S];  b.nutMass_t = 75e-6;   % autorotation nut
spR = buildSeedParams(b, cfgR);
results = {};

%% A. Added-mass rate (switched ON)
rng(3);  errA1 = 0;
for k = 1:200
    x = [zeros(3,1); randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100];
    x(4:7) = x(4:7)/norm(x(4:7));
    [~, im] = seed6DOFODE3D(0, x, spR);
    R  = quatToRotm(x(4:7));  mp = getMassProperties(0, spR);
    Ai = R*mp.A_trans*R.';    wI = R*x(11:13);   v = x(8:10);
    ref = cross(wI, Ai*v) - Ai*cross(wI, v);
    errA1 = max(errA1, norm(im.F_addedMassRate - ref)/max(norm(ref),eps));
end
results(end+1,:) = {'A1 implementation == formula', errA1 < 1e-12, sprintf('max rel err %.2e', errA1)};

sol = ode45(@(t,y) seed6DOFODE3D(t,y,spR), [0 6], [0;0;0;1;0;0;0;0;0;0;0;0;0], ...
            odeset('RelTol',1e-10,'AbsTol',1e-12));
mp  = getMassProperties(0, spR);  Ab = mp.A_trans;
Ai  = @(t) quatToRotm(deval(sol,t,4:7)/norm(deval(sol,t,4:7))) * Ab * ...
           quatToRotm(deval(sol,t,4:7)/norm(deval(sol,t,4:7))).';
h = 1e-6;  errA2 = 0;
for t = linspace(2, 5.9, 60)
    v  = deval(sol,t,8:10);
    fd = (Ai(t+h) - Ai(t-h))/(2*h) * v;
    [~, im] = seed6DOFODE3D(t, deval(sol,t), spR);
    errA2 = max(errA2, norm(fd - im.F_addedMassRate)/max(norm(fd),eps));
end
results(end+1,:) = {'A2 formula == finite difference of R A_b R''', errA2 < 1e-6, sprintf('max rel err %.2e', errA2)};

x = [zeros(3,1); axisAngleToQuat([1;1;0]/sqrt(2), 0.7); 0.3;-1.2;0.4; 5;30;-8];
spOff = spR;  spOff.enableAddedMassRate = false;
[dOn,  imOn]  = seed6DOFODE3D(0, x, spR);
[dOff, ~   ]  = seed6DOFODE3D(0, x, spOff);
R  = quatToRotm(x(4:7));  mp = getMassProperties(0, spR);
Meff = mp.M*eye(3) + R*mp.A_trans*R.';
predicted = -(Meff \ imOn.F_addedMassRate);
errA3 = norm((dOn(8:10)-dOff(8:10)) - predicted)/norm(predicted);
results(end+1,:) = {'A3 a_on - a_off == -(M+A_i)\(Adot v)', errA3 < 1e-10, sprintf('rel err %.2e', errA3)};

xNoSpin = x;  xNoSpin(11:13) = 0;
[~, im0]  = seed6DOFODE3D(0, xNoSpin, spR);
[~, imOf] = seed6DOFODE3D(0, x, spOff);
okA4 = norm(im0.F_addedMassRate) == 0 && norm(imOf.F_addedMassRate) == 0;
results(end+1,:) = {'A4 zero with omega = 0, and when switched off', okA4, ...
    sprintf('|F| no-spin %.1e, off %.1e', norm(im0.F_addedMassRate), norm(imOf.F_addedMassRate))};

% A5: Kirchhoff body-frame form (the APW form) vs this model's inertial form.
%   Body frame, M_b = M*I + A_b:   M_b vdot_b = F_b - w x (M_b v_b)
%   and the inertial acceleration is a = R (vdot_b + w x v_b).
% Same state, same net force F (aero + gravity) from the RHS intermediates.
rng(5);  errOn = 0;  errOff = 0;
for k = 1:100
    xk = [zeros(3,1); randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100];
    xk(4:7) = xk(4:7)/norm(xk(4:7));
    [dk, imk]   = seed6DOFODE3D(0, xk, spR);
    [dk0, ~]    = seed6DOFODE3D(0, xk, spOff);
    Rk = quatToRotm(xk(4:7));  mpk = getMassProperties(0, spR);
    Mb = mpk.M*eye(3) + mpk.A_trans;
    wb = xk(11:13);   vb = Rk.' * xk(8:10);
    Fb = Rk.' * imk.F_total_inertial;                  % aero + gravity, body frame
    aKirch = Rk * ( Mb \ (Fb - cross(wb, Mb*vb)) + cross(wb, vb) );
    errOn  = max(errOn,  norm(dk(8:10)  - aKirch)/norm(aKirch));
    errOff = max(errOff, norm(dk0(8:10) - aKirch)/norm(aKirch));
end
okA5 = errOn < 1e-12 && errOff > 1e-3;
results(end+1,:) = {'A5 ON == Kirchhoff/APW body-frame form; OFF is not', okA5, ...
    sprintf('ON rel err %.1e, OFF rel err %.1e', errOn, errOff)};

spDefault3 = buildSeedParams(b, cfg3);   spDefaultP = buildSeedParams(b, cfg);
okA6 = spDefault3.enableAddedMassRate && ~spDefaultP.enableAddedMassRate;
results(end+1,:) = {'A6 ON by default for shape3d, OFF for frozen planar', okA6, ''};

%% M. Added-mass (Munk) moment -- v x (A v), Kirchhoff's rotational partner
% Tested switched ON. This is APW eq. (6.3)'s (m11 - m22) vx' vy' term, which the
% model had no counterpart for at all.
cfgM = cfg3;  cfgM.enableAddedMassMoment = true;
spM  = buildSeedParams(b, cfgM);
spMoff = spM;  spMoff.enableAddedMassMoment = false;

rng(21);  errM1 = 0;
for k = 1:200
    xk = [zeros(3,1); randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100];
    xk(4:7) = xk(4:7)/norm(xk(4:7));
    [~, im] = seed6DOFODE3D(0, xk, spM);
    Rk  = quatToRotm(xk(4:7));   mpk = getMassProperties(0, spM);
    vb  = Rk.' * xk(8:10);
    ref = cross(vb, mpk.A_trans * vb);
    errM1 = max(errM1, norm(im.tau_addedMass - ref)/max(norm(ref),eps));
end
results(end+1,:) = {'M1 implementation == v x (A v)', errM1 < 1e-12, sprintf('max rel err %.2e', errM1)};

% M2: the 2D reduction that ties it to APW eq. (6.3). In-plane motion only
% (v_span = 0) must give a span-axis moment of exactly (A_chord - A_normal) vx vy,
% which is their (m11 - m22) vx' vy' with m11 <-> A_chord, m22 <-> A_normal.
mpM = getMassProperties(0, spM);   A = mpM.A_trans;
errM2 = 0;
for k = 1:100
    vb = [randn; randn; 0];                       % purely in-plane body velocity
    tauSolve = -cross(vb, A*vb);                  % what rigidBody6DOF adds to tau
    apw = (A(1,1) - A(2,2)) * vb(1) * vb(2);      % APW (6.3), first term
    errM2 = max(errM2, abs(tauSolve(3) - apw)/max(abs(apw),eps));
    errM2 = max(errM2, norm(tauSolve(1:2)));      % and nothing about the other axes
end
results(end+1,:) = {'M2 2D limit == APW eq. (6.3) (m11-m22)vx''vy''', errM2 < 1e-12, ...
    sprintf('max rel err %.2e; A_chord %.3e, A_normal %.3e kg', errM2, A(1,1), A(2,2))};

% M3: enters the rotational solve with the right sign.
xk = [zeros(3,1); axisAngleToQuat([1;-2;0.5]/norm([1;-2;0.5]), 1.1); 0.8;-1.7;0.6; 4;-12;7];
[dOnM,  imOnM] = seed6DOFODE3D(0, xk, spM);
[dOffM, ~    ] = seed6DOFODE3D(0, xk, spMoff);
mpk = getMassProperties(0, spM);   Ieff = mpk.I_G + mpk.A_rot;
predM = -(Ieff \ imOnM.tau_addedMass);
errM3 = norm((dOnM(11:13) - dOffM(11:13)) - predM)/norm(predM);
results(end+1,:) = {'M3 alpha_on - alpha_off == -(I+A_rot)\(v x A v)', errM3 < 1e-10, ...
    sprintf('rel err %.2e', errM3)};

% M4: vanishes along a principal axis, and when switched off.
okM4 = true;
for ax = 1:3
    vb = zeros(3,1);  vb(ax) = 2.5;
    okM4 = okM4 && norm(cross(vb, A*vb)) < 1e-18;
end
[~, imM0] = seed6DOFODE3D(0, xk, spMoff);
okM4 = okM4 && norm(imM0.tau_addedMass) == 0;
results(end+1,:) = {'M4 zero on a principal axis, and when switched off', okM4, ...
    sprintf('|tau| off = %.1e', norm(imM0.tau_addedMass))};

% M5: THE PHYSICAL CLAIM -- it destabilises an edge-on fall. Sliding along the span
% (v_z = V) with a small normal velocity +/-eps, the moment must GROW that
% misalignment, i.e. turn the plate broadside-on. Body-frame kinematics give
% dv_y/dt = omega_x v_z, so the growth condition is sign(tau_x) = sign(eps).
V = 3.0;   okM5 = true;   tauSample = 0;
for ep = [-0.3, -0.05, 0.05, 0.3]
    vb  = [0; ep; V];
    tau = -cross(vb, A*vb);                       % as added to the torque
    okM5 = okM5 && sign(tau(1)) == sign(ep) && abs(tau(2)) < 1e-18;
    if ep > 0; tauSample = tau(1); end
end
expectM5 = (A(2,2) - A(3,3)) * 0.3 * V;           % (A_normal - A_span) eps V
okM5 = okM5 && abs(tauSample - expectM5)/abs(expectM5) < 1e-12;
results(end+1,:) = {'M5 destabilises an edge-on fall (turns it broadside)', okM5, ...
    sprintf('tau_x = %+.3e N*m at eps = 0.3 m/s, V = %.1f m/s', tauSample, V)};

okM6 = spDefault3.enableAddedMassMoment && ~spDefaultP.enableAddedMassMoment;
results(end+1,:) = {'M6 ON by default for shape3d, OFF for frozen planar', okM6, ''};

%% B. Edge drag
A_edge = th*c;  Cd = 1.2;  zhat = [0;0;1];  arm = [0.001; 0; 0.02];
F  = @(v) computeEdgeDrag(v, arm, zhat, A_edge, Cd, rho);
f0 = F([0.7; -1.3; 0]);
results(end+1,:) = {'B1 zero for zero spanwise velocity', norm(f0) == 0, sprintf('|F| = %.1e', norm(f0))};
v1 = [0.2; -0.4; 0.9];
okB2 = norm(F(2*v1) - 4*F(v1)) < 1e-18 && norm(F(-v1) + F(v1)) < 1e-18;
results(end+1,:) = {'B2 quadratic and odd in v_s', okB2, ''};
okB3 = true;
for k = 1:200
    vr = randn(3,1);  fr = F(vr);
    okB3 = okB3 && (fr.'*zhat)*(vr.'*zhat) <= 0;
end
results(end+1,:) = {'B3 always opposes the spanwise velocity', okB3, '200 random states'};
fUnit = F([0;0;1]);  expect = 0.5*rho*Cd*A_edge;
okB4 = abs(norm(fUnit) - expect)/expect < 1e-14;
results(end+1,:) = {'B4 magnitude == 1/2 rho C_d (t c) v^2', okB4, ...
    sprintf('%.4e N at 1 m/s (expect %.4e)', norm(fUnit), expect)};
bc = b;  bc.nutPos_t = [0;0;0];  spc = buildSeedParams(bc, cfg3);
xs = [zeros(3,1); 1;0;0;0; 0;0;1.0; 0;0;0];
[~, ims] = seed6DOFODE3D(0, xs, spc);
stripZ = sum(ims.F_strip_body(3,:));
okB5 = abs(ims.F_aero_body(3) - ims.F_edge(3)) < 1e-18 && abs(stripZ) < 1e-18 && ims.F_edge(3) < 0;
results(end+1,:) = {'B5 pure spanwise slide: net F_z == edge drag', okB5, ...
    sprintf('strips F_z %.1e, edge F_z %.4e N', stripZ, ims.F_edge(3))};
spcOff = spc;  spcOff.enableEdgeDrag = false;
[~, imsOff] = seed6DOFODE3D(0, xs, spcOff);
okB6 = norm(imsOff.F_edge) == 0 && abs(imsOff.F_aero_body(3)) < 1e-18;
results(end+1,:) = {'B6 switched off -> zero, and the slide is unresisted', okB6, ...
    sprintf('net F_z off = %.1e', imsOff.F_aero_body(3))};

%% L. LEV vortex lift
% L1: planform constants
L = levPlanformConstants(S/c);
okL1a = abs(L.AR-3.3333)<1e-3 && abs(L.Kp-3.557)<1e-3 && abs(L.Ki-0.0955)<1e-4 && abs(L.Kv-2.349)<1e-3;
Lbig = levPlanformConstants(1e6);   Lsmall = levPlanformConstants(1e-3);
okL1b = abs(Lbig.Kp - 2*pi)/(2*pi) < 1e-5 && abs(Lsmall.Kp - pi*1e-3/2)/(pi*1e-3/2) < 1e-5;
KvScan = arrayfun(@(ar) getfield(levPlanformConstants(ar), 'Kv'), logspace(-1, 2, 200));
okL1c = all(KvScan > 0);
spL = buildSeedParams(b, cfg3);
okL1d = abs(spL.lev.Kv - L.Kv) < 1e-15;                  % builder uses total span / mean chord
results(end+1,:) = {'L1 constants: seed values, Helmbold limits, Kv>0, builder', ...
    okL1a && okL1b && okL1c && okL1d, ...
    sprintf('AR %.2f  Kp %.3f  Ki %.4f  Kv %.3f', L.AR, L.Kp, L.Ki, L.Kv)};

% L2: coefficient shape and sign, ungated
% "Ungated" = the GEOMETRIC definition at r = 0, which gives G = 1 exactly. (The
% kinematic gate has no ungated state: G -> 1 needs Omega -> inf.) The gate itself
% is L3 (geometric) and L3b/L3c (kinematic).
lev = spL.lev;   ch = c;   dzz = 0.005;
levG = lev;  levG.rossbyDefinition = 'geometric';
aGrid = deg2rad(-179.5:0.5:179.5);
okL2 = true;  errMag = 0;
for aa = aGrid
    [~, ~, inf1] = computeLEVForce(cos(aa), sin(aa), aa, ch, dzz, 0, 0, 0, levG, 'colocated', rho);
    co = computeAeroCoeffs(aa, []);
    errMag = max(errMag, abs(abs(inf1.CT_lev) - lev.Kv*sin(aa)^2*abs(cos(aa))));
    if abs(co.CT) > 1e-9 && abs(inf1.CT_lev) > 1e-9
        okL2 = okL2 && sign(inf1.CT_lev) == sign(co.CT);
    end
end
zeroAt = @(aa) computeLEVForce(cos(aa), sin(aa), aa, ch, dzz, 0, 0, 0, levG, 'colocated', rho);
[F0,~,i0] = zeroAt(0);  [F9,~,i9] = zeroAt(pi/2);  [Fm9,~,im9] = zeroAt(-pi/2);
okL2 = okL2 && errMag < 1e-14 && abs(i0.CT_lev) < 1e-15 && abs(i9.CT_lev) < 1e-15 ...
       && abs(im9.CT_lev) < 1e-15 && norm(F0) < 1e-18;
results(end+1,:) = {'L2 zero at 0/+-90 deg, Kv sin^2|cos|, APW sign', okL2, ...
    sprintf('max |mag err| %.1e over 719 angles', errMag)};

% L3: Rossby gate, GEOMETRIC definition (Ro = r/c)
Gof = @(r) getfield(nthOut3(@computeLEVForce, 1, 1, pi/4, ch, dzz, r, 0, 0, levG, 'colocated', rho), 'G');
rG  = linspace(0, 10*ch, 200);  Gs = arrayfun(Gof, rG);
okL3 = abs(Gof(0) - 1) < 1e-15 && abs(Gof(lev.RoCrit*ch) - 0.5) < 1e-15 && all(diff(Gs) < 0);
results(end+1,:) = {'L3 geometric gate: G(0)=1, G(Ro_crit)=1/2, decreasing', okL3, ...
    sprintf('G at Ro = 1, 3, 5: %.3f, %.3f, %.3f', Gof(ch), Gof(3*ch), Gof(5*ch))};

% L3b: Rossby gate, KINEMATIC definition (Ro = |v_ip|/(Omega c)) -- the default.
% The claim that justifies it: under PURE REVOLUTION, v_ip = Omega*r, so it must
% reproduce the geometric gate exactly. Plus the limits that make it well-posed for
% the integrator: no rotation -> shut; no flow -> open; neither -> finite, not NaN.
Gkin = @(vip, Om, r) nthOut3(@computeLEVForce, vip*cos(pi/4), vip*sin(pi/4), pi/4, ...
                             ch, dzz, r, Om, 0, lev, 'colocated', rho);
errRed = 0;  errRo = 0;
rng(11);
for k = 1:200
    rr = 0.5*ch + 8*ch*rand;   Om = 0.5 + 200*rand;      % pure revolution: v_ip = Om*r
    ik = Gkin(Om*rr, Om, rr);  ig = getfield(nthOut3(@computeLEVForce, 1, 1, pi/4, ch, dzz, ...
                                    rr, 0, 0, levG, 'colocated', rho), 'G');
    errRed = max(errRed, abs(ik.G - ig)/ig);
    errRo  = max(errRo,  abs(ik.Ro - rr/ch)/(rr/ch));
end
iNoSpin = Gkin(2.0, 0,   0.02);       % translating, not revolving
iNoFlow = Gkin(0,   40,  0.02);       % revolving, momentarily no in-plane flow
iNeither = Gkin(0,  0,   0.02);       % both zero: 0/0 must not appear
iHalf   = Gkin(lev.RoCrit*40*ch, 40, 0.02);                  % Ro == Ro_crit
vScan   = linspace(0.01, 12, 300);
Gscan   = arrayfun(@(vv) getfield(Gkin(vv, 25, 0.02), 'G'), vScan);
okL3b = errRed < 1e-12 && errRo < 1e-12 ...
        && iNoSpin.G == 0 && isinf(iNoSpin.Ro) ...
        && abs(iNoFlow.G - 1) < 1e-15 && iNoFlow.Ro == 0 ...
        && iNeither.G == 0 && isfinite(iNeither.G) ...
        && abs(iHalf.G - 0.5) < 1e-14 && all(diff(Gscan) < 0);
results(end+1,:) = {'L3b kinematic gate == geometric under pure revolution; limits', okL3b, ...
    sprintf('reduction rel err %.1e (Ro %.1e); Omega=0 -> G %.0f, v=0 -> G %.0f', ...
            errRed, errRo, iNoSpin.G, iNoFlow.G)};

% L3c: the same thing INSIDE the RHS, on the case the geometric gate got wrong.
%   (i) a seed translating with NO rotation gets zero LEV under 'kinematic' and a
%       large one under 'geometric' -- this is the defect the switch fixes.
%  (ii) a TWISTED seed in pure revolution (v = 0, omega about body y) has
%       v_ip = Omega*r on every strip, so the two definitions must agree there.
%       Twist is what makes alpha (and so the force) nonzero in that state.
spKin = spL;  spKin.enableLEV = true;                        % builder default: kinematic
spGeo = spKin;  spGeo.lev.rossbyDefinition = 'geometric';
xTrans = [zeros(3,1); 1;0;0;0; 0.4;-2.0;0.1; 0;0;0];         % falling, not spinning
[~, iKt] = seed6DOFODE3D(0, xTrans, spKin);
[~, iGt] = seed6DOFODE3D(0, xTrans, spGeo);
bTw = b;  bTw.nutPos_t = [0;0;0];  bTw.twist = @(z) (deg2rad(20)/(S/2))*z;
spTk = buildSeedParams(bTw, cfg3);  spTk.enableLEV = true;
spTg = spTk;  spTg.lev.rossbyDefinition = 'geometric';
xRev = [zeros(3,1); 1;0;0;0; 0;0;0; 0;45;0];                 % pure revolution about body y
[~, iKr] = seed6DOFODE3D(0, xRev, spTk);
[~, iGr] = seed6DOFODE3D(0, xRev, spTg);
redRHS = max(abs(iKr.G_lev - iGr.G_lev) ./ iGr.G_lev);
% (This seed's nut is at 1.2*S, so its CoM is outboard and the geometric gate is
%  partly closed on the far strips -- Ro = r/c reaches ~3.4 there. The point is that
%  it is open AT ALL on a seed that is not revolving, not that it is open fully.)
okL3c = all(iKt.G_lev == 0) && norm(iKt.F_lev_body(:)) == 0 ...
        && all(iGt.G_lev > 0) && max(iGt.G_lev) > 0.9 && norm(iGt.F_lev_body(:)) > 0 ...
        && redRHS < 1e-12 && norm(iKr.F_lev_body(:)) > 0;
results(end+1,:) = {'L3c in RHS: no rotation -> no LEV; pure revolution -> == geometric', okL3c, ...
    sprintf('translating G: kin %.3f vs geo %.3f-%.3f; revolving gate rel err %.1e, |F_lev| %.1e N', ...
            max(iKt.G_lev), min(iGt.G_lev), max(iGt.G_lev), redRHS, norm(iKr.F_lev_body(:)))};

% L4: it is lift -- perpendicular to the in-plane velocity
okL4 = true;
for k = 1:200
    vc = randn;  vn = randn;  aa = atan2(vn, vc);
    Fl = computeLEVForce(vc, vn, aa, ch, dzz, 0.01, 0, 0, levG, 'colocated', rho);
    okL4 = okL4 && abs(Fl(1)*vc + Fl(2)*vn) < 1e-15*max(1,norm(Fl)) && Fl(3) == 0;
end
results(end+1,:) = {'L4 force is lift: F . v_inplane = 0', okL4, '200 random states'};

% L5: application points
aa = deg2rad(35);  co = computeAeroCoeffs(aa, []);
[~, xC] = computeLEVForce(cos(aa), sin(aa), aa, ch, dzz, 0, 0, co.l_cp_frac, levG, 'colocated', rho);
[~, xF] = computeLEVForce(cos(aa), sin(aa), aa, ch, dzz, 0, 0, co.l_cp_frac, levG, 'forward',   rho);
aR = deg2rad(145);  % reversed flow: the -chord edge leads
[~, xFr] = computeLEVForce(cos(aR), sin(aR), aR, ch, dzz, 0, 0, 0, levG, 'forward', rho);
[~, xF9] = computeLEVForce(0, 1, pi/2, ch, dzz, 0, 0, 0, levG, 'forward', rho);
okL5 = abs(xC - co.l_cp_frac*ch) < 1e-18 && abs(xF - lev.lambdaV*(ch/2)*cos(aa)) < 1e-18 ...
       && xF > 0 && xFr < 0 && abs(xF9) < 1e-18;
threw = false;
try; computeLEVForce(1, 0, 0, ch, dzz, 0, 0, 0, levG, 'fwd', rho); catch; threw = true; end
results(end+1,:) = {'L5 colocated = l_cp*c; forward upwind, 0 at 90 deg', okL5 && threw, ...
    sprintf('x/c at 35 deg: colocated %+.3f, forward %+.3f; at 145 deg forward %+.3f', ...
            xC/ch, xF/ch, xFr/ch)};

% L6: wiring inside the RHS
spLon  = spL;  spLon.enableLEV = true;   spLon.levApplicationPoint = 'colocated';
spLfwd = spLon;                          spLfwd.levApplicationPoint = 'forward';
spLnofield = rmfield(spL, {'enableLEV','levApplicationPoint','lev'});
rng(9);  okL6 = true;  errF = 0;  errT = 0;  dT = 0;
for k = 1:100
    xk = [zeros(3,1); randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*60];
    xk(4:7) = xk(4:7)/norm(xk(4:7));
    [dOffL, iOffL] = seed6DOFODE3D(0, xk, spL);
    dNone          = seed6DOFODE3D(0, xk, spLnofield);
    okL6 = okL6 && isequal(dOffL, dNone) && all(iOffL.CT_lev == 0);
    [~, iOn] = seed6DOFODE3D(0, xk, spLon);
    [~, iFw] = seed6DOFODE3D(0, xk, spLfwd);
    mpk = getMassProperties(0, spL);
    tauLev = zeros(3,1);
    for i = 1:numel(spL.strips.chord)
        cD = spL.strips.chordDir(:,i);
        pG = [spL.strips.xgc_body(i); spL.strips.ygc_body(i); spL.strips.zgc_body(i)];
        tauLev = tauLev + cross(pG + iOn.xApp_lev(i)*cD - mpk.c, iOn.F_lev_body(:,i));
    end
    % relative to the totals: (a + b) - a recovers b only to rounding of a
    errF = max(errF, norm((iOn.F_aero_body - iOffL.F_aero_body) - sum(iOn.F_lev_body,2)) ...
                     / norm(iOn.F_aero_body));
    errT = max(errT, norm((iOn.tau_body - iOffL.tau_body) - tauLev) / norm(iOn.tau_body));
    okL6 = okL6 && isequal(iOn.F_aero_body, iFw.F_aero_body);   % same force exactly
    dT   = max(dT, norm(iOn.tau_body - iFw.tau_body));
end
okL6 = okL6 && errF < 1e-12 && errT < 1e-12 && dT > 0;
results(end+1,:) = {'L6 in RHS: off bit-identical; on adds exactly F_lev, r x F_lev', okL6, ...
    sprintf('rel force err %.1e, rel torque err %.1e; colocated vs forward torques differ by up to %.1e N*m', ...
            errF, errT, dT)};

% L7: defaults and overrides
okL7 = ~spL.enableLEV && strcmp(spL.levApplicationPoint, 'colocated') && spL.lev.RoCrit == 3 ...
       && spL.lev.p == 4 && spL.lev.lambdaV == 0.5 ...
       && strcmp(spL.lev.rossbyDefinition, 'kinematic');
cK = cfg3;  cK.lev = struct('rossbyDefinition', 'geometric');   % the toggle is settable
sK = buildSeedParams(b, cK);
okL7 = okL7 && strcmp(sK.lev.rossbyDefinition, 'geometric') && sK.lev.Kv == spL.lev.Kv;
threwRo = false;
try; cR = cfg3; cR.lev = struct('rossbyDefinition','kinetic'); buildSeedParams(b, cR); catch; threwRo = true; end
okL7 = okL7 && threwRo;
cO = cfg3;  cO.lev = struct('Kp', 3.0);                        % Kv must follow eq. (4)
sO = buildSeedParams(b, cO);
okL7 = okL7 && abs(sO.lev.Kv - (3.0 - 9.0*sO.lev.Ki)) < 1e-15;
cA = cfg3;  cA.lev = struct('AR', 4.38, 'RoCrit', 2.5);         % AR re-derives Kp, Ki, Kv
sA = buildSeedParams(b, cA);   LA = levPlanformConstants(4.38);
okL7 = okL7 && abs(sA.lev.Kv - LA.Kv) < 1e-15 && sA.lev.RoCrit == 2.5;
threw = false;
try; cB = cfg3; cB.levApplicationPoint = 'front'; buildSeedParams(b, cB); catch; threw = true; end
results(end+1,:) = {'L7 defaults off/colocated/kinematic; cfg.lev keeps Kv = Kp - Kp^2 Ki', okL7 && threw, ...
    sprintf('Kp=3 -> Kv %.3f;  AR 4.38 -> Kv %.3f;  bad names rejected', sO.lev.Kv, sA.lev.Kv)};

%% C. Planar untouched
spP = buildSeedParams(b, cfg);
okC = ~spP.enableAddedMassRate && ~isfield(spP,'enableEdgeDrag') && ~isfield(spP,'enableLEV') ...
      && ~isfield(spP,'lev');
results(end+1,:) = {'C  planar: rate OFF, no edge-drag or LEV switch', okC, ''};

%% Report
fprintf('\n');
allOk = true;
for k = 1:size(results,1)
    fprintf('  %-58s %-5s %s\n', results{k,1}, pf(results{k,2}), results{k,3});
    allOk = allOk && results{k,2};
end
if allOk
    fprintf('\nPASS: added-mass rate, edge drag and LEV behave as derived, inside and outside the RHS.\n');
else
    fprintf(2, '\nFAIL: see the checks above.\n');
end

function s = pf(ok);  if ok; s = 'PASS'; else; s = 'FAIL'; end;  end
function out = nthOut3(fn, varargin)
% Return the THIRD output of fn (computeLEVForce's info struct).
    [~, ~, out] = fn(varargin{:});
end
