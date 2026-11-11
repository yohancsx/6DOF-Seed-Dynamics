function res = testRealSeedLoader()
% TESTREALSEEDLOADER  The measured-seed import (loadRealSeeds / realSeedToBsp) and
%   the builder's explicit-geometry input (bsp.stripGeom in setupSeedShape3D).
%
% CASES
%   R1  Contract: every seed in seedMassModel.mat builds (buildSeedParams), passes
%       validateSeedParams, the builder passes stripGeom through UNCHANGED, and the
%       RHS is finite at random 3D states.
%   R2  Geometry: sum(chord.*dz) = the plates' fitted area; sum(dz) = their length;
%       cutting each plate into 3 sub-strips keeps the area and its centroid.
%   R3  Mass, calibration 'both': M = measured mass, and the CoM and inertia are the
%       file's own calibrated values (out.calibrated), moved to body axes by P.
%   R4  Mass, calibration 'nutOnly': M = measured, wing = its predicted mass, nut =
%       the rest, inertia symmetric positive definite; and with predicted = calibrated
%       part masses it reduces EXACTLY to 'both'.
%   R5  Flat limit: 12 synthetic flat rectangular plates + a point nut at the origin,
%       through realSeedToBsp and stripGeom, against the polyshape builder's flat
%       rectangle (12 strips, same nut): strips, mass properties and the RHS.
%   R6  Resolution: 1 vs 4 sub-strips per plate on three seeds; the aerodynamic
%       force at a fixed state is reported (WARN above 5%).
% If the data file is absent (it lives outside the repo), the test is skipped with
% a WARN, so the suite still runs on machines without it.
%
% OUTPUT
%   res : struct array (.name, .pass, .detail).

    here = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(here,'seedPaths.m')) && ~strcmp(here,fileparts(here)); here = fileparts(here); end
    addpath(here);  Pth = seedPaths(here);
    res = struct('name', {}, 'pass', {}, 'detail', {});
    if ~isfile(Pth.seedMassModel)
        res = item('real-seed loader', true, ['WARN skipped: ' Pth.seedMassModel ' not found']);
        report(res);  return
    end
    P   = [0 1 0; 0 0 1; 1 0 0];
    cfg = struct('rhoFluid', 1.225, 'g', 9.81, 'enableAddedMass3D', true);
    D   = load(Pth.seedMassModel, 'out', 'TP');

    % ============================== R1 / R2 / R3 =================================
    seeds = loadRealSeeds(Pth.seedMassModel);
    s3    = loadRealSeeds(Pth.seedMassModel, struct('subStrips', 3));
    rng(31);  okC = true;  passThru = true;  finite = true;  eA = 0;  eL = 0;  eSub = 0;
    eM = 0;  eCom = 0;  eI = 0;  msg = '';
    for k = 1:numel(seeds)
        o  = D.out([D.out.seedNum] == seeds(k).seedNum);
        try
            sp = buildSeedParams(seeds(k).bsp, cfg);   validateSeedParams(sp);
        catch err
            okC = false;  msg = sprintf('%s: %s', seeds(k).name, err.message);  continue
        end
        sg = seeds(k).bsp.stripGeom;
        passThru = passThru && isequaln(sp.strips, sg.strips) && isequaln(sp.massParams, sg.massParams);
        for j = 1:5
            q = randn(4,1);  q = q/norm(q);
            xd = seed6DOFODE3D(0, [zeros(3,1); q; 2*randn(3,1); 20*randn(3,1)], sp);
            finite = finite && all(isfinite(xd));
        end
        st = sp.strips;   pl = o.plates.plates;
        eA = max(eA, abs(sum(st.chord.*st.dz) / (sum([pl.area])*1e-6) - 1));
        eL = max(eL, abs(sum(st.dz) / (sum([pl.length])*1e-3) - 1));
        t3 = s3(k).bsp.stripGeom.strips;
        c1 = [st.xgc_body; st.ygc_body; st.zgc_body] * st.area.' / sum(st.area);
        c3 = [t3.xgc_body; t3.ygc_body; t3.zgc_body] * t3.area.' / sum(t3.area);
        eSub = max([eSub, abs(sum(t3.area)/sum(st.area) - 1), norm(c3 - c1)/0.001]);
        C = o.calibrated;  mp = sp.massParams;
        eM   = max(eM,   abs(mp.M_total / (o.measuredMass*1e-3) - 1));
        eCom = max(eCom, norm(mp.com_t - P*C.com(:)*1e-3) / 1e-3);
        eI   = max(eI,   norm(mp.I_G_t - P*C.I*P.'*1e-9) / norm(C.I*1e-9));
    end
    res(end+1) = item(sprintf('R1 contract: %d seeds build and run', numel(seeds)), ...
                      okC && passThru && finite, sprintf('validate %d, stripGeom passed through %d, RHS finite %d %s', ...
                      okC, passThru, finite, msg));
    res(end+1) = item('R2 geometry: area, length, sub-strips', max([eA eL eSub]) < 1e-12, ...
                      sprintf('area rel %.1e, length rel %.1e, 1 vs 3 sub-strips (area, centroid/mm) %.1e', eA, eL, eSub));
    res(end+1) = item('R3 mass (both) = file''s calibrated values', max([eM eCom eI]) < 1e-12, ...
                      sprintf('M rel %.1e, CoM %.1e mm, I_G rel %.1e (body axes)', eM, eCom, eI));

    % ============================== R4 ===========================================
    sN = loadRealSeeds(Pth.seedMassModel, struct('calibration', 'nutOnly'));
    eM = 0;  eW = 0;  pd = true;  eRed = 0;
    for k = 1:numel(sN)
        o  = D.out([D.out.seedNum] == sN(k).seedNum);
        parts = D.TP(D.TP.seed == string(o.name), :);
        mp = sN(k).bsp.stripGeom.massParams;  me = sN(k).meta;
        eM = max(eM, abs(mp.M_total/(o.measuredMass*1e-3) - 1));
        eW = max(eW, abs(me.wingMass / (sum(parts.predMass_g(parts.body=="plate"))*1e-3) - 1));
        Ig = mp.I_G_t;   pd = pd && norm(Ig - Ig.') < 1e-15*norm(Ig) && all(eig((Ig+Ig.')/2) > 0);
        % predicted = calibrated part masses -> nutOnly must equal both
        p2 = parts;  p2.predMass_g = p2.mass_g;   o2 = o;  o2.measuredMass = sum(p2.mass_g);
        bB = realSeedToBsp(o2, p2, struct('calibration', 'both', 'subStrips', 1));
        bN = realSeedToBsp(o2, p2, struct('calibration', 'nutOnly', 'subStrips', 1));
        mB = bB.stripGeom.massParams;  mN = bN.stripGeom.massParams;
        eRed = max([eRed, abs(mN.M_total/mB.M_total - 1), norm(mN.com_t - mB.com_t)/1e-3, ...
                    norm(mN.I_G_t - mB.I_G_t)/norm(mB.I_G_t)]);
    end
    nutFr = arrayfun(@(s) s.meta.nutFrac, sN);   nutFrB = arrayfun(@(s) s.meta.nutFrac, seeds);
    res(end+1) = item('R4 mass (nutOnly)', max([eM eW]) < 1e-12 && pd && eRed < 1e-12, sprintf( ...
        'M rel %.1e, wing = predicted %.1e, PD %d, reduces to both %.1e;  nut fraction %.2f-%.2f (both: %.2f-%.2f)', ...
        eM, eW, pd, eRed, min(nutFr), max(nutFr), min(nutFrB), max(nutFrB)));

    % ============================== R5 ===========================================
    S = 0.050;  c = 0.015;  th = 0.002;  rhoA = 65*th;  nP = 12;  mNut = 75e-6;
    b.seedShape = polyshape([-S/2, S/2, S/2, -S/2], [-c/2, -c/2, c/2, c/2]);
    b.seedDensity = rhoA;  b.seedThickness = th;  b.numStrips = nP;  b.tSamples = 0;
    b.nutPos_t = [0;0;0];  b.nutMass_t = mNut;
    ref = buildSeedParams(b, cfg);
    [o, parts] = flatPlateSeed(S*1e3, c*1e3, th*1e3, rhoA, nP, mNut);
    bE = realSeedToBsp(o, parts, struct('calibration', 'both', 'subStrips', 1));
    bE.seedThickness = th;
    syn = buildSeedParams(bE, cfg);
    a = ref.strips;  e = syn.strips;
    dS = max([abs(a.chord - e.chord), abs(a.dz - e.dz), abs(a.xgc_body - e.xgc_body), ...
              abs(a.zgc_body - e.zgc_body), abs(a.ygc_body - e.ygc_body), ...
              max(abs(a.chordDir - e.chordDir)), max(abs(a.normalDir - e.normalDir)), ...
              max(abs(a.spanDir - e.spanDir))]);
    dMc = norm(ref.massParams.com_t - syn.massParams.com_t);
    dMI = norm(ref.massParams.I_G_t - syn.massParams.I_G_t) / norm(ref.massParams.I_G_t);
    dMM = abs(ref.massParams.M_total - syn.massParams.M_total) / ref.massParams.M_total;
    rng(32);  dX = 0;
    for j = 1:20
        q = randn(4,1);  q = q/norm(q);  x = [zeros(3,1); q; 2*randn(3,1); 20*randn(3,1)];
        xa = seed6DOFODE3D(0, x, ref);  xe = seed6DOFODE3D(0, x, syn);
        dX = max(dX, max(abs(xa - xe) ./ max(abs(xa), 1e-3)));
    end
    res(end+1) = item('R5 flat limit = polyshape builder', max([dS/1e-3, dMc/1e-3, dMI, dMM, dX]) < 1e-9, ...
        sprintf('strips %.1e m, CoM %.1e m, I_G rel %.1e, M rel %.1e, RHS rel %.1e', dS, dMc, dMI, dMM, dX));

    % ============================== R6 ===========================================
    pick = [5 24 29];  dF = zeros(size(pick));
    for k = 1:numel(pick)
        a1 = buildSeedParams(loadRealSeeds(Pth.seedMassModel, struct('seeds', pick(k))).bsp, cfg);
        a4 = buildSeedParams(loadRealSeeds(Pth.seedMassModel, struct('seeds', pick(k), 'subStrips', 4)).bsp, cfg);
        x  = [zeros(3,1); axisAngleToQuat([0;0;1], deg2rad(30)); 0.3; -1.2; 0.1; 2; 5; 30];
        [~, i1] = seed6DOFODE3D(0, x, a1);   [~, i4] = seed6DOFODE3D(0, x, a4);
        dF(k) = norm(i4.F_aero_body - i1.F_aero_body) / norm(i1.F_aero_body);
    end
    w = '';  if max(dF) > 0.05; w = 'WARN '; end
    res(end+1) = item('R6 resolution: 1 vs 4 sub-strips', true, sprintf( ...
        '%saero force change at a tumbling state: seeds %s -> %s', w, mat2str(pick), mat2str(100*dF, 2) + "%"));

    report(res);
end

% =========================================================================
function [o, parts] = flatPlateSeed(Smm, cmm, tmm, rhoA, nP, mNut)
% A flat rectangle as 12 'measured' plates in the file's own layout (mesh axes,
% mm / g / g mm^2): x = span, y = chord, z = normal. Thin-lamina inertias.
    xs = linspace(-Smm/2, Smm/2, nP + 1);
    pl = struct('e1', {}, 'e2', {}, 'n', {}, 'corners', {}, 'length', {}, 'area', {}, ...
                'twistDeg', {}, 'dihedralDeg', {}, 'thickness', {});
    rows = {};
    for i = 1:nP
        L = xs(i+1) - xs(i);
        pl(i) = struct('e1', [1 0 0], 'e2', [0 1 0], 'n', [0 0 1], ...
            'corners', [xs(i) -cmm/2 0; xs(i) cmm/2 0; xs(i+1) cmm/2 0; xs(i+1) -cmm/2 0], ...
            'length', L, 'area', L*cmm, 'twistDeg', 0, 'dihedralDeg', 0, 'thickness', tmm);
        m = rhoA * (L*cmm*1e-6) * 1e3;                         % g
        rows(end+1, :) = {"plate", m, m, (xs(i)+xs(i+1))/2, 0, 0, ...
            m*cmm^2/12, m*L^2/12, m*(cmm^2 + L^2)/12, 0, 0, 0}; %#ok<AGROW>
    end
    rows(end+1, :) = {"nut", mNut*1e3, mNut*1e3, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    parts = cell2table(rows, 'VariableNames', {'body', 'mass_g', 'predMass_g', ...
        'centroidLumped_x', 'centroidLumped_y', 'centroidLumped_z', 'Ixx_lumped_gmm2', ...
        'Iyy_lumped_gmm2', 'Izz_lumped_gmm2', 'Ixy_lumped_gmm2', 'Ixz_lumped_gmm2', 'Iyz_lumped_gmm2'});
    parts.body = string(parts.body);
    o = struct('name', 'flatSynthetic', 'seedNum', 0, 'confidence', "synthetic", ...
               'measuredMass', sum(parts.mass_g), 'calibration', 1);
    o.plates.plates = pl;
end
function report(res)
    for k = 1:numel(res)
        if ~res(k).pass; tag = 'FAIL'; elseif startsWith(res(k).detail, 'WARN'); tag = 'WARN'; else; tag = 'PASS'; end
        fprintf('  %-44s %s  %s\n', res(k).name, tag, res(k).detail);
    end
end
function r = item(name, pass, detail)
    r = struct('name', name, 'pass', logical(pass), 'detail', char(detail));
end
