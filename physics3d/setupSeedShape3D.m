function seedParamsFull = setupSeedShape3D(seedParamsIn)
% SETUPSEEDSHAPE3D  Build a NON-PLANAR (3D-shape) seed's geometry + mass.
%
% The 'shape3d' analogue of setupSeedShapeAndMass. It produces the same strip +
% mass-property fields the ODE consumes, PLUS the per-strip 3D contract fields
% (ygc_body and the local chord/normal/span frames) that let seed6DOFODE3D place
% and orient each strip out of the body x-z plane (twist / camber / dihedral).
%
% STATUS: supports flat plates, TWIST (spanwise geometric pitch about the span
% axis), and spanwise CURVATURE (dihedral bending about the chord axis, strips
% leaving y=0). It wraps setupSeedShapeAndMass for the planar strips, then sets
% each strip's out-of-plane position + local frame from the requested twist and
% curvature and retags 'shape3d'. With neither, the frames are identity and the
% seed is a strict superset of the planar seed (the flat-equivalence regression).
%
% MASS / INERTIA: twist keeps every strip centre at y=0 and only re-orients a
% thin plate about its own mid-chord, so the flat-plate mass properties stay
% exact and are used unchanged. CURVATURE moves mass off the x-z plane, so the
% CoM and inertia are recomputed from the 3D strip positions + orientations
% (see curvedMassProperties below): the wing centroid picks up a y component and
% each strip's local inertia is rotated into its own frame before the
% parallel-axis sum. Added mass is still the flat-plate form (a documented
% approximation, gentle-curvature only, pending its own generalisation).
%
% TWIST INPUT: baseSeedParams.twist = scalar (uniform pitch, rad) / length-M
% vector / handle @(z) (pitch vs spanwise position) -- e.g. anti-symmetric
% @(z) k*z. Absent -> no twist.
% CURVATURE INPUT: baseSeedParams.curvature = the dihedral tangent angle phi vs
% ARC LENGTH s, same forms -- e.g. @(s) k*s for a circular-arc bowl, or a scalar
% for a uniform tilt. Absent -> flat (no bending). Twist and curvature compose.
%
% INPUT
%   seedParamsIn : struct with .baseSeedParams (a flat polyshape planform + nut
%                  mass/position over time, plus optional .twist / .curvature).
% OUTPUT
%   seedParamsFull : the planar build, plus
%                    .strips.ygc_body  (1xM, out-of-plane strip position, body y)
%                    .strips.chordDir/normalDir/spanDir (3xM strip frame, body)
%                    .strips.twist     (1xM applied pitch per strip, rad)
%                    .strips.dihedral  (1xM applied tangent angle per strip, rad)
%                    .model = 'shape3d'
%                  (.strips.zgc_body/z_body are remapped to the physical spanwise
%                  position on the curve when curvature is set.)

    bsp = seedParamsIn.baseSeedParams;

    % --- Planar base: identical strips + mass to the flat-plate builder -----
    seedParamsFull = setupSeedShapeAndMass(seedParamsIn);

    M    = numel(seedParamsFull.strips.chord);
    sArc = seedParamsFull.strips.zgc_body;       % flat span coordinate = ARC LENGTH
                                                  % along the (bent) profile (1xM)

    % --- TWIST: per-strip geometric pitch theta(s) about the LOCAL SPAN axis ---
    % Twist rotates each strip's frame about the span axis by theta, changing its
    % incidence while KEEPING its centre in place. An anti-symmetric theta(s)
    % turns some of each strip's lift into a chordwise force at a spanwise arm ->
    % a centred-CoM seed spins up (the samara autorotation mechanism).
    theta = resolveTwist(bsp, sArc, M);          % 1xM pitch per strip (rad); 0 -> flat

    % --- CURVATURE: spanwise dihedral tangent angle phi(s) about the CHORD axis -
    % The flat span coordinate sArc is treated as ARC LENGTH; a curvature bends the
    % profile in the body (y,z) plane. phi(s) is the local tangent (dihedral) angle;
    % integrating dz=cos(phi)ds, dy=sin(phi)ds gives each strip's PHYSICAL (z,y)
    % centre on the curve, anchored so s=0 -> (0,0). Strip widths (dz) are the
    % arc-length segments, so they are UNCHANGED by bending. Absent -> flat
    % (phi=0, y=0, positions unchanged), so twist-only and planar stay bit-identical.
    %
    hasCurvature = isfield(bsp, 'curvature') && ~isempty(bsp.curvature);
    if hasCurvature
        [phi, zCurve, yCurve] = curveFromProfile(bsp, sArc, M);
        seedParamsFull.strips.zgc_body = zCurve;    % physical spanwise position (body z)
        seedParamsFull.strips.z_body   = zCurve;
        ygc = yCurve;
    else
        phi = zeros(1, M);
        ygc = zeros(1, M);
    end

    % --- Per-strip frame = dihedral(phi, about chord) o twist(theta, about span).
    % Columns of R_strip = Rd(phi)*Rz(theta), in body coords. Reduces to the twist
    % frame when phi=0 and to identity when theta=phi=0 (so flat/twist unchanged).
    ct = cos(theta);   st = sin(theta);   cp = cos(phi);   sp = sin(phi);   % 1xM
    seedParamsFull.strips.chordDir  = [ct;         cp.*st;   -sp.*st];
    seedParamsFull.strips.normalDir = [-st;        cp.*ct;   -sp.*ct];
    seedParamsFull.strips.spanDir   = [zeros(1,M); sp;        cp];
    seedParamsFull.strips.ygc_body  = ygc;
    seedParamsFull.strips.twist     = theta;     % record profiles (rad)
    seedParamsFull.strips.dihedral  = phi;

    % --- MASS / INERTIA ----------------------------------------------------
    % Twist keeps every strip centre at y = 0 and only re-orients a thin plate
    % about its own mid-chord, so the flat-plate mass properties from
    % setupSeedShapeAndMass remain exact -- leave them alone (this is also what
    % keeps flat/twist bit-identical to the planar model). CURVATURE genuinely
    % moves mass off the x-z plane, so recompute the CoM and inertia from the 3D
    % strip positions + orientations.
    if hasCurvature
        seedParamsFull.massParams = curvedMassProperties(seedParamsFull, bsp);

        % --- Span-force hack OFF by default for a CURVED seed --------------
        % computeSpanForce is planar heritage: it fakes the body-z force that a
        % FLAT seed's strips can never produce (their normals are all body-y).
        % Once the strips are tilted out of plane they generate that spanwise
        % force from geometry, so leaving the hack on DOUBLE-COUNTS it. Twist
        % alone does not tilt the normals out of the x-y plane, so it keeps the
        % planar default. Override via cfg.enableSpanForce if you want to A/B it.
        seedParamsFull.enableSpanForce = false;
    end

    % --- Drop the switches the 3D RHS no longer honours --------------------
    % setupSeedShapeAndMass stamps the full planar switch set. seed6DOFODE3D no
    % longer reads these four -- the terms they gated (the span torque with its
    % migrating CoP and reduced-frequency attenuation, and the Tx / Ty whole-seed
    % spin-damping torques) were removed by the Sep-2026 audit. Leaving them on
    % the struct would advertise knobs that silently do nothing, so strip them:
    % a shape3d seedParams carrying one of these names is a stale config, and
    % now looks like one. The planar model still has all of them.
    deadSwitches = {'enableSpanTorque', 'enableSpanCOPMigration', ...
                    'enableSpanTorqueAttenuation', 'enableTxDamping', ...
                    'enableNormalSpinDamping'};
    for iDead = 1:numel(deadSwitches)
        if isfield(seedParamsFull, deadSwitches{iDead})
            seedParamsFull = rmfield(seedParamsFull, deadSwitches{iDead});
        end
    end

    % --- Retag the model (setupSeedShapeAndMass stamped it 'planar') --------
    seedParamsFull.model = 'shape3d';
end


% =========================================================================
% LOCAL: 3D mass properties for a seed whose strips leave the body x-z plane
% =========================================================================
function mp = curvedMassProperties(sp, bsp)
% Lumps each strip as a thin uniform rectangle (mass m_i = rho_A*A_i, chord
% extent c_i, span-width w_i) sitting at its 3D centroid r_i with orientation
% R_i = [chordDir normalDir spanDir], and adds the nut as a point mass:
%
%   CoM      r_G = ( sum_i m_i r_i + m_n r_n ) / M
%   local    I_i' = (1/12) m_i * diag( w_i^2, c_i^2 + w_i^2, c_i^2 )
%                   [strip axes: chord, normal, span; I_22' by the
%                    perpendicular-axis theorem for a lamina]
%   rotate   I_i^loc = R_i * I_i' * R_i'
%   Steiner  I_G = sum_i [ I_i^loc + m_i( |d_i|^2 I3 - d_i d_i' ) ]
%                    + m_n( |d_n|^2 I3 - d_n d_n' ),   d = r - r_G
%
% With R_i = I and y_i = 0 this collapses term-for-term onto the planar
% builder's flat-plate computation. I_G_dot uses the same finite-difference
% scheme as the planar builder (only r_G and the nut term move in time; the
% shape itself is fixed in the body frame).

    s    = sp.strips;
    rhoA = bsp.seedDensity;              % areal density (kg/m^2)
    tS   = bsp.tSamples;                 % time samples (as supplied)
    numT = numel(tS);
    M    = numel(s.chord);

    % --- Per-strip mass, centroid, and body-frame local inertia ------------
    m_i = rhoA * s.area(:).';                        % 1xM (actual clipped areas)
    r_i = [s.xgc_body; s.ygc_body; s.zgc_body];      % 3xM strip centroids
    c   = s.chord(:).';                              % chord extent per strip
    w   = s.dz(:).';                                 % span-width (arc length) per strip

    Iloc = zeros(3, 3, M);
    for i = 1 : M
        Ri = [s.chordDir(:,i), s.normalDir(:,i), s.spanDir(:,i)];   % strip -> body
        Ip = (1/12) * m_i(i) * diag([ w(i)^2, c(i)^2 + w(i)^2, c(i)^2 ]);
        Iloc(:,:,i) = Ri * Ip * Ri.';
    end

    m_w    = sum(m_i);                   % wing mass
    S_wing = r_i * m_i(:);               % 3x1 first moment of the wing, sum_i m_i r_i

    % --- CoM + inertia at each mass time sample (the nut may move) ---------
    com_t = zeros(3, numT);
    I_G_t = zeros(3, 3, numT);
    for k = 1 : numT
        m_n = bsp.nutMass_t(k);
        r_n = bsp.nutPos_t(:, k);
        rG  = (S_wing + m_n * r_n) / (m_w + m_n);
        com_t(:, k) = rG;

        Ik = zeros(3);
        for i = 1 : M
            d  = r_i(:, i) - rG;
            Ik = Ik + Iloc(:,:,i) + m_i(i) * ((d.' * d) * eye(3) - d * d.');
        end
        dn = r_n - rG;
        Ik = Ik + m_n * ((dn.' * dn) * eye(3) - dn * dn.');
        I_G_t(:, :, k) = Ik;
    end

    % --- d/dt of I_G: central differences inside, one-sided at the ends ----
    I_G_dot_t = zeros(3, 3, numT);
    if numT > 1
        for k = 1 : numT
            if k == 1
                I_G_dot_t(:,:,k) = (I_G_t(:,:,2) - I_G_t(:,:,1)) / (tS(2) - tS(1));
            elseif k == numT
                I_G_dot_t(:,:,k) = (I_G_t(:,:,end) - I_G_t(:,:,end-1)) / (tS(end) - tS(end-1));
            else
                I_G_dot_t(:,:,k) = (I_G_t(:,:,k+1) - I_G_t(:,:,k-1)) / (tS(k+1) - tS(k-1));
            end
        end
    end

    mp = struct('tSamples', tS, 'com_t', com_t, 'I_G_t', I_G_t, ...
                'I_G_dot_t', I_G_dot_t, 'M_total', m_w + bsp.nutMass_t(1));
end


% =========================================================================
% LOCAL: curvature profile -> per-strip tangent angle + physical (z,y) centres
% =========================================================================
function [phi, zNew, yNew] = curveFromProfile(bsp, sCenters, M)
% bsp.curvature (same forms as bsp.twist) gives the dihedral tangent angle phi as
% a function of ARC LENGTH s: a scalar (uniform tilt -> a shallow V), a length-M
% vector (per strip), or a handle @(s) (e.g. @(s) k*s for a circular-arc bowl).
% Returns phi at the strip centres and the physical (z,y) centre of each strip on
% the integrated curve (arc-length s=0 anchored to the origin).
    cv = bsp.curvature;
    ev = @(sq) evalProfile(cv, sq, sCenters, M);   % phi at query arc-lengths
    phi = ev(sCenters);

    sGrid = linspace(min(sCenters), max(sCenters), 1001);
    phiG  = ev(sGrid);
    zG = cumtrapz(sGrid, cos(phiG));               % dz = cos(phi) ds
    yG = cumtrapz(sGrid, sin(phiG));               % dy = sin(phi) ds
    zG = zG - interp1(sGrid, zG, 0, 'linear', 'extrap');   % anchor s=0 -> (0,0)
    yG = yG - interp1(sGrid, yG, 0, 'linear', 'extrap');
    zNew = interp1(sGrid, zG, sCenters, 'linear', 'extrap');
    yNew = interp1(sGrid, yG, sCenters, 'linear', 'extrap');
end

function v = evalProfile(cv, sq, sCenters, M)
% Evaluate a scalar / length-M vector / handle profile at arc-lengths sq.
    if isa(cv, 'function_handle')
        v = arrayfun(@(s) cv(s), sq(:).');
    elseif isscalar(cv)
        v = double(cv) * ones(1, numel(sq));
    elseif numel(cv) == M
        v = interp1(sCenters, double(cv(:)).', sq, 'linear', 'extrap');
    else
        error('setupSeedShape3D:badCurvature', ...
              'bsp.curvature must be a scalar, a length-%d vector, or a handle @(s).', M);
    end
end


% =========================================================================
% LOCAL: resolve the twist input into a 1xM per-strip pitch (rad)
% =========================================================================
function theta = resolveTwist(bsp, zgc, M)
% bsp.twist may be: absent/empty (flat, theta = 0); a scalar (uniform pitch on
% every strip); a length-M vector (per-strip pitch, strip order); or a function
% handle @(z) giving the pitch at spanwise position z (m), evaluated at each
% strip's spanwise centre -- e.g. a linear anti-symmetric twist
% @(z) k*z/(spanLength/2).
    if ~isfield(bsp, 'twist') || isempty(bsp.twist)
        theta = zeros(1, M);
        return
    end
    tw = bsp.twist;
    if isa(tw, 'function_handle')
        theta = arrayfun(@(z) tw(z), zgc(:).');          % 1xM
    elseif isscalar(tw)
        theta = double(tw) * ones(1, M);
    elseif numel(tw) == M
        theta = double(tw(:)).';
    else
        error('setupSeedShape3D:badTwist', ...
              'bsp.twist must be a scalar, a length-%d vector, or a function handle @(z).', M);
    end
end
