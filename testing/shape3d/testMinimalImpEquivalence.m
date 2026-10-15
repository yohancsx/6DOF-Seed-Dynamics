function res = testMinimalImpEquivalence()
% TESTMINIMALIMPEQUIVALENCE  The 3D model, confined to 2D, against Olivia's
%   published 2D reference model (Olivia Code/minimal_imp.m).
%
% minimal_imp is the validated Andersen-Pesavento-Wang / Pomerenk-Ristroph plate
% model that every sectional coefficient here comes from. This checks that the 3D
% strip model, flown in the plane, IS that model -- not a re-derivation of it.
%
% THE REFERENCE IS CALLED, NOT COPIED. Its Model function is extracted VERBATIM
% from minimal_imp.m at run time into a temporary folder and called there. Olivia's
% files are never edited, and no copy of her code lives in this repo.
%
% MATCHING THE TWO MODELS
%   units        minimal_imp is CGS (cm, g, s); ours is SI. Accelerations compared
%                in CGS.
%   per span     minimal_imp is per unit depth. Ours is a uniform rectangle of
%                span S cut into strips; forces and masses all scale with S, so
%                accelerations compare directly and forces are divided by S.
%   fluid        rhoFluid = 1000 kg/m^3 (her p = 1 g/cm^3, water); g = 9.80.
%   mass         massParams overridden (test only) to her m*S, I*S and CoM offset
%                lcm along the chord: M_total, I_G(3,3) and com_t(1).
%   buoyancy     ours has none, so mp = 0 in her call (her Tb and (m-mp) terms).
%   added mass   ours matches her m1 = 0, m2 = pi rho (c/2)^2, Ia = pi rho c^4/128 +
%                m2 lcm^2 by construction (getAddedMass). Her omega*v cross terms
%                are our added-mass RATE term (ON; testNewPhysicsTerms A5). She has
%                NO Munk moment (m11-m22)vx'vy', so enableAddedMassMoment = false.
%   other        LEV off. Edge drag at its default: the spanwise velocity is zero
%                in the plane, so it is inert.
%   state        her [x y theta vx' vy' omega] <-> our r, q = rotation theta about
%                world Z, v = R [vx'; vy'; 0], omega = [0; 0; omega]. Her vx', vy'
%                are the CoM velocity (her alpha and v_cg are formed at the
%                geometric centre, vy' - omega*lcm), exactly as our strips are.
%
% CASES
%   E1  lcm = 0, 500 random states: (d vx'/dt, d vy'/dt, d omega/dt) agree.
%       Every CoM-offset term (Tb, Tcr, -m2 omega^2 lcm, the explicit d omega/dt
%       in her F_yp) vanishes, so the two must be the same equations.
%   E2  lcm ~= 0: the aerodynamic force and torque agree exactly. Her aero force is
%       recovered from her own outputs by removing her (known) inertial terms; her
%       torque from d omega/dt (her rotation equation is torque/(I+Ia), as ours is
%       for planar motion).
%   E3  lcm ~= 0: the full translational equations. Here the two models are
%       STRUCTURALLY different -- she adds -m2 omega^2 lcm and m2 lcm d omega/dt
%       (the added mass of an off-centre plate seen from its CoM), ours has no
%       translation-rotation added-mass coupling. Reported as WARN, never failed.
%   E4  Trajectories from her own initial condition (theta = -30 deg, 5 cm/s,
%       5 s): at lcm = 0 they must coincide until chaos separates them (the
%       agreement horizon is reported); at lcm = 0.1 cm the descent rate and the
%       rotation are compared, and differences reported (they include E3's).
%   E5  Planarity: a 2D release stays in the plane (z, v_z, omega_x, omega_y = 0).
%
% OUTPUT
%   res : struct array (.name, .pass, .detail).

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  P = seedPaths(here);

    % --- her model, extracted verbatim -----------------------------------------
    src = fileread(fullfile(P.root, 'Olivia Code', 'minimal_imp.m'));
    k0  = regexp(src, 'function\s+y\s*=\s*Model\s*\(', 'once');
    assert(~isempty(k0), 'minimal_imp.m: Model function not found');
    tmp = tempname;  mkdir(tmp);
    fid = fopen(fullfile(tmp, 'Model.m'), 'w');  fwrite(fid, src(k0:end));  fclose(fid);
    addpath(tmp);
    cleanup = onCleanup(@() rmTemp(tmp)); %#ok<NASGU>

    % --- her parameters (CGS), from minimal_imp.m ----------------------------------
    l = 2.54;  m = 1;  I = 0.5;  p = 1;  g = 980;
    m2 = pi*p*l^2/4;
    % --- our seed: a uniform rectangle in SI -------------------------------------
    S = 0.10;  c = l/100;
    res = struct('name', {}, 'pass', {}, 'detail', {});

    % ============================== E1 =========================================
    sp0 = herSeed(0);
    rng(11);  N = 500;  err = zeros(N, 3);
    for j = 1:N
        st = [0; 0; 2*pi*rand - pi; 40*randn; 40*randn; 15*randn];
        [a_ours, ~] = ours(sp0, st, 0);
        a_hers = Model(st, 0, I, l, m, 0, p, g);
        err(j,:) = abs(a_ours - a_hers(4:6)).' ./ [g g max(abs(a_hers(6)), 1)];
    end
    res(end+1) = item('E1 lcm = 0: equations of motion', max(err(:)) < 1e-10, ...
        sprintf('max rel err  dvx %.1e  dvy %.1e  domega %.1e  (500 states)', max(err)));

    % ============================== E2 / E3 ====================================
    lcm = 0.1;   spL = herSeed(lcm);   Ia = pi*p*l^4*(1 + 32*lcm^2/l^2)/128;
    rng(12);  eF = zeros(N,2);  eT = zeros(N,1);  eA = zeros(N,2);
    for j = 1:N
        st = [0; 0; 2*pi*rand - pi; 40*randn; 40*randn; 15*randn];
        [a_ours, F_ours] = ours(spL, st, lcm);
        y  = Model(st, lcm, I, l, m, 0, p, g);
        th = st(3);  vx = st(4);  vy = st(5);  w = st(6);  T = (I + Ia)*y(6);
        % her aero force: her EOM with her inertial and gravity terms removed
        Fx = m*y(4) - (m + m2)*w*vy + m2*w^2*lcm + m*g*sin(th);
        Fy = (m + m2)*y(5) + m*w*vx - m2*T/(I + Ia)*lcm + m*g*cos(th);
        eF(j,:) = abs(F_ours - [Fx Fy]) ./ max(norm([Fx Fy]), m*g);
        eT(j)   = abs(a_ours(3) - y(6)) / max(abs(y(6)), 1);
        eA(j,:) = abs(a_ours(1:2) - y(4:5)).' / g;
    end
    res(end+1) = item('E2 lcm = 0.1 cm: aero force and torque', max([eF(:); eT]) < 1e-10, ...
        sprintf('max rel err  F %.1e  torque (via domega) %.1e', max(eF(:)), max(eT)));
    res(end+1) = item('E3 lcm = 0.1 cm: full translational EOM', true, sprintf( ...
        ['%sdvx max %.2e g, median %.2e g;  dvy max %.2e g, median %.2e g  ' ...
         '(her -m2 w^2 lcm and m2 lcm dw/dt terms; ours has no such coupling)'], ...
        warnIf(max(eA(:)) > 1e-10), max(eA(:,1)), median(eA(:,1)), max(eA(:,2)), median(eA(:,2))));

    % ============================== E4 / E5 ====================================
    x0h = [0; 0; -pi/6; 5; 0; 0];   tq = linspace(0, 5, 2001).';
    tight = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    for lc = [0 0.1]
        [~, Xh] = ode45(@(t, x) Model(x, lc, I, l, m, 0, p, g), tq, x0h, tight);
        sp = herSeed(lc);
        R0 = rotz2(x0h(3));
        x0 = [zeros(3,1); axisAngleToQuat([0;0;1], x0h(3)); R0*[x0h(4); x0h(5); 0]/100; 0; 0; x0h(6)];
        [~, Xo] = ode45(@(t, x) seed6DOFODE3D(t, x, sp), tq, x0, tight);
        th_o = unwrap(2*atan2(Xo(:,7), Xo(:,4)));      % rotation about Z from the quaternion
        dPos = hypot(100*Xo(:,1) - Xh(:,1), 100*Xo(:,2) - Xh(:,2));   % cm
        dTh  = abs(th_o - Xh(:,3));
        bad  = find(dPos > 0.01*l | dTh > 0.01, 1);      % 1% of a chord, or 0.6 deg
        if isempty(bad); horizon = tq(end); else; horizon = tq(bad); end
        vd_h = -(Xh(end,2) - Xh(1,2)) / tq(end);   vd_o = -100*(Xo(end,2) - Xo(1,2)) / tq(end);
        turns_h = (Xh(end,3) - Xh(1,3)) / (2*pi);   turns_o = (th_o(end) - th_o(1)) / (2*pi);
        det = sprintf(['agree to 1%% chord for %.2f of 5 s;  descent %.2f vs her %.2f cm/s;  ' ...
                       'rotation %+.2f vs her %+.2f turns'], horizon, vd_o, vd_h, turns_o, turns_h);
        if lc == 0
            res(end+1) = item('E4 lcm = 0: trajectory', horizon >= 1.0, det); %#ok<AGROW>
        else
            res(end+1) = item('E4 lcm = 0.1 cm: trajectory', true, ...
                              [warnIf(abs(vd_o - vd_h) > 0.05*abs(vd_h)) det]); %#ok<AGROW>
        end
        plan = max([abs(Xo(:,3)); abs(Xo(:,10)); abs(Xo(:,11)); abs(Xo(:,12))]);
        res(end+1) = item(sprintf('E5 lcm = %.1f cm: stays planar', lc), plan < 1e-9, ...
                          sprintf('max |z|, |v_z|, |w_x|, |w_y| = %.1e', plan)); %#ok<AGROW>
    end

    % --- report ------------------------------------------------------------------
    for k = 1:numel(res)
        if ~res(k).pass; tag = 'FAIL'; elseif startsWith(res(k).detail, 'WARN'); tag = 'WARN'; else; tag = 'PASS'; end
        fprintf('  %-40s %s  %s\n', res(k).name, tag, res(k).detail);
    end
    if all([res.pass])
        fprintf('PASS: the 3D model flown in the plane reproduces minimal_imp (E3 and E4b report known differences).\n');
    else
        fprintf(2, 'FAIL: the 3D model in 2D departs from minimal_imp (see above).\n');
    end

    % =========================================================================
    function sp = herSeed(lcmCm)
    % Her plate as our seed: uniform rectangle, her mass, inertia and CoM offset.
        b.seedShape     = polyshape([-S/2, S/2, S/2, -S/2], [-c/2, -c/2, c/2, c/2]);
        b.seedDensity   = 1;   b.seedThickness = 1e-4;   b.numStrips = 10;
        b.tSamples = 0;   b.nutPos_t = [0;0;0];   b.nutMass_t = 0;
        cfg = struct('rhoFluid', 1000, 'g', 9.80, 'shapeModel', 'shape3d', ...
                     'enableAddedMassRate', true, 'enableAddedMassMoment', false, ...
                     'enableLEV', false);
        sp = buildSeedParams(b, cfg);
        mp = sp.massParams;
        mp.M_total      = (m/1e3*100) * S;                  % g/cm -> kg/m, times span
        mp.com_t(:,1)   = [lcmCm/100; 0; 0];
        Ig              = mp.I_G_t(:,:,1);
        Ig(3,3)         = (I*1e-3*1e-4/1e-2) * S;           % g cm^2/cm -> kg m^2/m
        Ig(1,3) = 0;  Ig(3,1) = 0;  Ig(2,3) = 0;  Ig(3,2) = 0;
        mp.I_G_t(:,:,1) = Ig;
        mp.I_G_dot_t    = zeros(size(mp.I_G_dot_t));
        sp.massParams   = mp;
    end
    function [aCgs, Fcgs] = ours(sp, st, ~)
    % Our model at her state: body-frame accelerations (CGS) and aero force per span.
        R  = rotz2(st(3));
        u  = [st(4); st(5); 0] / 100;                       % body velocity, m/s
        w  = [0; 0; st(6)];
        x  = [zeros(3,1); axisAngleToQuat([0;0;1], st(3)); R*u; w];
        [xd, im] = seed6DOFODE3D(0, x, sp);
        du = R.'*xd(8:10) - cross(w, u);                    % d/dt of the body-frame velocity
        aCgs = [100*du(1); 100*du(2); xd(13)];
        Fcgs = (im.F_aero_body(1:2).' / S) * 1e3;           % N/m -> dyn/cm
    end
end

function R = rotz2(th)
    R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
end
function s = warnIf(tf)
    if tf; s = 'WARN '; else; s = ''; end
end
function r = item(name, pass, detail)
    r = struct('name', name, 'pass', logical(pass), 'detail', detail);
end
function rmTemp(d)
    rmpath(d);  rmdir(d, 's');
end
