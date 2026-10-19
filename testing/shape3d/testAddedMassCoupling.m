function res = testAddedMassCoupling()
% TESTADDEDMASSCOUPLING  The 2026-10-16 added-mass terms: the off-centre coupling
%   block A_coup (getAddedMass), Kirchhoff's 6x6 solve (rigidBody6DOF), and the
%   'residual' Munk mode. Their 2D limit against minimal_imp is checked in
%   testMinimalImpEquivalence (E3, E4, E6); this covers the 3D side.
%
% CASES
%   K1  A_coup from the flat-plate branch equals the per-strip-normal branch on a
%       flat seed, and its pitch entry is the 2D -m2*lcm for a chordwise offset.
%   K2  With the CoM at the added-mass centroid (A_coup = 0) the coupled 6x6 path
%       reproduces the uncoupled path (to round-off): it adds nothing else.
%   K3  Ideal fluid (no aero, no gravity), Hou plate with an off-centre nut, 2 s of
%       3D tumbling: the symmetric form conserves the kinetic energy
%       1/2 [V;W]' M6 [V;W]; the oneSided (minimal_imp) form does not.
%   K4  'residual' removes exactly the in-section (pitch) Munk moment and keeps the
%       spanwise (roll) part: at W = 0 on a flat plate, residual = (m2 v w, 0, 0).
%
% OUTPUT
%   res : struct array (.name, .pass, .detail).

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  seedPaths(here);
    res = struct('name', {}, 'pass', {}, 'detail', {});

    [pcfg, pbsp, toNutPos] = paperSeedConfig();
    mk = @(nut, cf) buildSeedParams(setfield(setfield(pbsp, 'nutPos_t', nut), ...
                                    'nutMass_t', pcfg.nutMass), cf); %#ok<SFLD>

    % ============================== K1 =========================================
    sp = mk(toNutPos(0.39, 0.17), pcfg);   mp = getMassProperties(0, sp);
    s3 = sp;  s3.enableAddedMass3D = true;
    [~, ~, Ac0] = getAddedMass(sp, mp.c);   [~, ~, Ac3] = getAddedMass(s3, mp.c);
    d13 = max(abs(Ac0(:) - Ac3(:))) / max(abs(Ac0(:)));
    m2t = sum(pi*sp.rhoFluid*(sp.strips.chord/2).^2 .* sp.strips.dz);
    xc  = mp.c(1) - sum(sp.strips.xgc_body .* sp.strips.chord.^2 .* sp.strips.dz) / ...
                    sum(sp.strips.chord.^2 .* sp.strips.dz);
    dPitch = abs(Ac0(2,3) + m2t*xc) / (m2t*abs(xc));
    res(end+1) = item('K1 A_coup: branches agree, 2D limit', d13 < 1e-12 && dPitch < 1e-12, ...
        sprintf('flat vs per-strip %.1e;  A_coup(2,3) vs -m2*lcm %.1e;  A_coup(2,:) = [%.2e %.0f %.2e] kg m', ...
                d13, dPitch, Ac0(2,1), Ac0(2,2), Ac0(2,3)));

    % ============================== K2 =========================================
    spC = mk([0;0;0], pcfg);                                 % nut at the plate centre
    [~, ~, AcC] = getAddedMass(spC, getMassProperties(0, spC).c);
    spU = spC;  spU.enableAddedMassCoupling = false;
    rng(21);  e = zeros(100,1);
    for j = 1:100
        q = randn(4,1);  q = q/norm(q);
        x = [zeros(3,1); q; 0.6*randn(3,1); 25*randn(3,1)];
        a = seed6DOFODE3D(0, x, spC);   b = seed6DOFODE3D(0, x, spU);
        e(j) = max(abs(a(8:13) - b(8:13)) ./ max(abs(b(8:13)), 1));
    end
    res(end+1) = item('K2 A_coup = 0: coupled path = uncoupled', max(e) < 1e-9, ...
        sprintf('max|A_coup| %.1e;  max rel diff %.1e over 100 3D states', max(abs(AcC(:))), max(e)));

    % ============================== K3 =========================================
    x0 = [zeros(3,1); [1;0.2;-0.3;0.1]/norm([1;0.2;-0.3;0.1]); 0.4*[1;-0.5;0.3]; [12;-20;30]];
    tight = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);
    drift = zeros(1,2);  forms = {'symmetric', 'oneSided'};
    for k = 1:2
        s = sp;  s.addedMassCouplingForm = forms{k};  s.addedMassMomentMode = 'full';
        mpk = getMassProperties(0, s);   mpk.tau_inSectionMunk = [0;0;0];
        [~, X] = ode45(@(t, xx) rigidBody6DOF(xx, mpk, [0;0;0], [0;0;0], 0), [0 2], x0, tight);
        M6 = [mpk.M*eye(3) + mpk.A_trans, mpk.A_coup; mpk.A_coup.', mpk.I_G + mpk.A_rot];
        E  = zeros(size(X,1),1);
        for j = 1:size(X,1)
            R  = quatToRotm(X(j,4:7).'/norm(X(j,4:7)));
            nu = [R.'*X(j,8:10).'; X(j,11:13).'];   E(j) = 0.5*nu.'*M6*nu;
        end
        drift(k) = max(abs(E/E(1) - 1));
    end
    res(end+1) = item('K3 ideal fluid: energy conserved', drift(1) < 1e-8 && drift(2) > 1e-4, ...
        sprintf('max |E/E0-1| over 2 s: symmetric %.1e;  oneSided %.1e (expected to drift)', drift));

    % ============================== K4 =========================================
    sR = sp;  sF = sp;  sF.addedMassMomentMode = 'full';
    rng(22);  e = zeros(50,1);
    for j = 1:50
        q = randn(4,1);  q = q/norm(q);  R = quatToRotm(q);
        x = [zeros(3,1); q; 0.6*randn(3,1); zeros(3,1)];
        [~, iR] = seed6DOFODE3D(0, x, sR);   [~, iF] = seed6DOFODE3D(0, x, sF);
        V = R.'*x(8:10);   mpj = getMassProperties(0, sR);   m2 = mpj.A_trans(2,2);
        expectF = [-m2*V(2)*V(3); 0; m2*V(1)*V(2)];     % V x (A V), the full term
        expectR = [-m2*V(2)*V(3); 0; 0];                % its spanwise (roll) part
        sc = max(norm(expectF), 1e-12);
        e(j) = max([norm(iF.tau_addedMass - expectF), norm(iR.tau_addedMass - expectR)]) / sc;
    end
    res(end+1) = item('K4 residual Munk keeps only the roll part', max(e) < 1e-10, ...
        sprintf('max rel err %.1e over 50 states (W = 0)', max(e)));

    % --- report ------------------------------------------------------------------
    for k = 1:numel(res)
        if res(k).pass; tag = 'PASS'; else; tag = 'FAIL'; end
        fprintf('  %-42s %s  %s\n', res(k).name, tag, res(k).detail);
    end
end

function r = item(name, pass, detail)
    r = struct('name', name, 'pass', logical(pass), 'detail', detail);
end
