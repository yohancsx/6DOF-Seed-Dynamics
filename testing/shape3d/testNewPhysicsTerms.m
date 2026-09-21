%% Regression: the phase-4 physics terms (added-mass rate + edge drag)
% Locks the two terms added to the shape3d model in phase 4, both ON by default.
% Run after any change to rigidBody6DOF, computeEdgeDrag or seed6DOFODE3D.
%
%   A. ADDED-MASS RATE (rigidBody6DOF, enableAddedMassRate)
%      A1  the implementation equals the formula  w x (A_i v) - A_i (w x v)
%      A2  the formula equals a finite difference of R*A_b*R' along a REAL
%          ode45 trajectory -- the physics, not just the algebra
%      A3  it enters the solve with the right sign:
%          a_on - a_off = -(M*I + A_i) \ (Adot_i v)
%      A4  it vanishes identically with no rotation, and when switched off
%   B. EDGE DRAG (computeEdgeDrag, enableEdgeDrag)
%      B1  zero for zero spanwise velocity
%      B2  quadratic:  F(2v) = 4 F(v);  odd: F(-v) = -F(v)
%      B3  always opposes the spanwise velocity
%      B4  magnitude matches 1/2 rho C_d (t c) v^2
%      B5  INSIDE THE RHS: for a pure spanwise slide the strips produce exactly
%          zero spanwise force, so the net body-z aero force IS the edge drag
%      B6  zero when switched off
%   C. PLANAR IS UNTOUCHED: the planar builder leaves the rate term OFF and has
%      no edge-drag switch at all.

root = 'C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics';
addpath(fullfile(root,'physics'), fullfile(root,'physics','helpers'), ...
        fullfile(root,'physics','aero'), fullfile(root,'physics','mass'), ...
        fullfile(root,'physics3d'), fullfile(root,'testing','helpers'));

S = 0.050;  c = 0.015;  th = 0.002;  rho = 1.225;
cfg  = struct('rhoFluid',rho,'g',9.81);
cfg3 = cfg;  cfg3.shapeModel = 'shape3d';
b.seedShape = polyshape([-S/2,S/2,S/2,-S/2],[-c/2,-c/2,c/2,c/2]);
b.seedDensity = 65*th;  b.seedThickness = th;  b.numStrips = 10;
b.tSamples = 0;  b.nutPos_t = [c;0;1.2*S];  b.nutMass_t = 75e-6;   % autorotation nut
sp3 = buildSeedParams(b, cfg3);
results = {};

%% A. Added-mass rate
% A1: implementation vs formula, on random states
rng(3);  errA1 = 0;
for k = 1:200
    x = [zeros(3,1); randn(4,1); (rand(3,1)-0.5)*6; (rand(3,1)-0.5)*100];
    x(4:7) = x(4:7)/norm(x(4:7));
    [~, im] = seed6DOFODE3D(0, x, sp3);
    R  = quatToRotm(x(4:7));  mp = getMassProperties(0, sp3);
    Ai = R*mp.A_trans*R.';    wI = R*x(11:13);   v = x(8:10);
    ref = cross(wI, Ai*v) - Ai*cross(wI, v);
    errA1 = max(errA1, norm(im.F_addedMassRate - ref)/max(norm(ref),eps));
end
results(end+1,:) = {'A1 implementation == formula', errA1 < 1e-12, sprintf('max rel err %.2e', errA1)};

% A2: formula vs finite difference along a real trajectory (term ON)
sol = ode45(@(t,y) seed6DOFODE3D(t,y,sp3), [0 6], [0;0;0;1;0;0;0;0;0;0;0;0;0], ...
            odeset('RelTol',1e-10,'AbsTol',1e-12));
mp  = getMassProperties(0, sp3);  Ab = mp.A_trans;
Ai  = @(t) quatToRotm(deval(sol,t,4:7)/norm(deval(sol,t,4:7))) * Ab * ...
           quatToRotm(deval(sol,t,4:7)/norm(deval(sol,t,4:7))).';
h = 1e-6;  errA2 = 0;
for t = linspace(2, 5.9, 60)
    v  = deval(sol,t,8:10);
    fd = (Ai(t+h) - Ai(t-h))/(2*h) * v;
    [~, im] = seed6DOFODE3D(t, deval(sol,t), sp3);
    errA2 = max(errA2, norm(fd - im.F_addedMassRate)/max(norm(fd),eps));
end
results(end+1,:) = {'A2 formula == finite difference of R A_b R''', errA2 < 1e-6, sprintf('max rel err %.2e', errA2)};

% A3: sign/entry into the translational solve
x = [zeros(3,1); axisAngleToQuat([1;1;0]/sqrt(2), 0.7); 0.3;-1.2;0.4; 5;30;-8];
spOff = sp3;  spOff.enableAddedMassRate = false;
[dOn,  imOn]  = seed6DOFODE3D(0, x, sp3);
[dOff, ~   ]  = seed6DOFODE3D(0, x, spOff);
R  = quatToRotm(x(4:7));  mp = getMassProperties(0, sp3);
Meff = mp.M*eye(3) + R*mp.A_trans*R.';
predicted = -(Meff \ imOn.F_addedMassRate);
errA3 = norm((dOn(8:10)-dOff(8:10)) - predicted)/norm(predicted);
results(end+1,:) = {'A3 a_on - a_off == -(M+A_i)\(Adot v)', errA3 < 1e-10, sprintf('rel err %.2e', errA3)};

% A4: vanishes with no rotation, and when off
xNoSpin = x;  xNoSpin(11:13) = 0;
[~, im0]  = seed6DOFODE3D(0, xNoSpin, sp3);
[~, imOf] = seed6DOFODE3D(0, x, spOff);
okA4 = norm(im0.F_addedMassRate) == 0 && norm(imOf.F_addedMassRate) == 0;
results(end+1,:) = {'A4 zero with omega = 0, and when switched off', okA4, ...
    sprintf('|F| no-spin %.1e, off %.1e', norm(im0.F_addedMassRate), norm(imOf.F_addedMassRate))};

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

fUnit = F([0;0;1]);
expect = 0.5*rho*Cd*A_edge*1^2;
okB4 = abs(norm(fUnit) - expect)/expect < 1e-14;
results(end+1,:) = {'B4 magnitude == 1/2 rho C_d (t c) v^2', okB4, ...
    sprintf('%.4e N at 1 m/s (expect %.4e)', norm(fUnit), expect)};

% B5: inside the RHS, a pure spanwise slide. Flat seed, NUT CENTRED so the area
% centroid sits on the CoM (no omega x r needed); omega = 0; level attitude so body
% z = world z. The strips see only spanwise flow and so produce exactly zero force.
bc = b;  bc.nutPos_t = [0;0;0];  spc = buildSeedParams(bc, cfg3);
xs = [zeros(3,1); 1;0;0;0; 0;0;1.0; 0;0;0];           % v = +1 m/s along world/body z
[~, ims] = seed6DOFODE3D(0, xs, spc);
stripZ = sum(ims.F_strip_body(3,:));
errB5 = abs(ims.F_aero_body(3) - ims.F_edge(3));
okB5 = errB5 < 1e-18 && abs(stripZ) < 1e-18 && ims.F_edge(3) < 0;
results(end+1,:) = {'B5 pure spanwise slide: net F_z == edge drag', okB5, ...
    sprintf('strips F_z %.1e, edge F_z %.4e N', stripZ, ims.F_edge(3))};

spcOff = spc;  spcOff.enableEdgeDrag = false;
[~, imsOff] = seed6DOFODE3D(0, xs, spcOff);
okB6 = norm(imsOff.F_edge) == 0 && abs(imsOff.F_aero_body(3)) < 1e-18;
results(end+1,:) = {'B6 switched off -> zero, and the slide is unresisted', okB6, ...
    sprintf('net F_z off = %.1e  (the unphysical zero this term removes)', imsOff.F_aero_body(3))};

%% C. Planar untouched
spP = buildSeedParams(b, cfg);
okC = ~spP.enableAddedMassRate && ~isfield(spP,'enableEdgeDrag');
results(end+1,:) = {'C  planar: rate OFF, no edge-drag switch', okC, ''};

%% Report
fprintf('\n');
allOk = true;
for k = 1:size(results,1)
    fprintf('  %-50s %-5s %s\n', results{k,1}, pf(results{k,2}), results{k,3});
    allOk = allOk && results{k,2};
end
if allOk
    fprintf('\nPASS: added-mass rate and edge drag behave as derived, inside and outside the RHS.\n');
else
    fprintf(2, '\nFAIL: see the checks above.\n');
end

function s = pf(ok);  if ok; s = 'PASS'; else; s = 'FAIL'; end;  end
