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
    end
    % (A curved seed used to switch the planar span-force hack off here, because
    % its tilted strips produce that force from geometry. The span force is now
    % retired from this model entirely -- see below -- so the special case is gone.
    % Its replacement, edge drag, is tip form drag: a different mechanism from the
    % strips' normal forces, so it does not double-count on a curved seed.)

    % --- Wing thickness: needed for the edge-drag frontal area ---------------
    % setupSeedShapeAndMass reads it but does not keep it; record it here rather
    % than change the frozen planar builder.
    seedParamsFull.seedThickness = bsp.seedThickness;

    % --- Physics switches for the shape3d model -------------------------------
    % Override any of these through cfg (buildSeedParams applies explicit overrides).
    %   enableEdgeDrag      ON  -- tip crossflow drag for spanwise sliding
    %                              (computeEdgeDrag). Replaces the span force and
    %                              removes a known-unphysical exact zero.
    %   enableAddedMassRate   ON -- the Adot*v term in the translational EOM
    %                              (rigidBody6DOF). The inertial-frame form of the
    %                              anisotropic added-mass terms the APW reference
    %                              model carries in its body-frame equations,
    %                              eqs. (6.1)-(6.2) (testNewPhysicsTerms check A5).
    %   enableAddedMassMoment ON -- v x (A*v), Kirchhoff's ROTATIONAL partner of
    %                              that term and APW eq. (6.3)'s (m11-m22)vx'vy'.
    %                              The moment that turns a body broadside-on to its
    %                              own motion, so it is what makes an edge-on fall
    %                              unstable. Also in rigidBody6DOF.
    %
    % Those two are ON by the rule that an added-mass term belongs in the model by
    % default when it is derived from first principles AND present in the reference
    % model. Both are: they are the two halves of Kirchhoff's pair (Lamb,
    % Hydrodynamics ch. VI), and both appear in APW 2005 -- the rate term as the
    % omega*v cross terms of eqs. (6.1)-(6.2), the moment as the first term of
    % eq. (6.3). Taking one without the other was taking half of Kirchhoff. They
    % also measurably improve the comparison with Hou et al. (2025): descent speed
    % moves from 1.2-2.4 to 1.0-1.4 times sqrt(sigma g/rho), where the published
    % prefactor is O(1). The FROZEN PLANAR builder keeps both OFF, so the planar
    % reference stays byte-identical to every earlier run.
    %   enableLEV           OFF -- leading-edge-vortex vortex lift (computeLEVForce).
    %                              Empirical, and not yet validated for this seed.
    %   levApplicationPoint 'colocated' -- where that vortex lift acts: at the APW
    %                              centre of pressure (Snyder & Lamar 1972), or
    %                              'forward' toward the leading edge (uncited).
    % (The LEV's other modelling switch, which Rossby number the gate reads, lives
    %  on the lev struct as .rossbyDefinition -- documented with the constants below.)
    seedParamsFull.enableEdgeDrag        = true;
    seedParamsFull.enableAddedMassRate   = true;
    seedParamsFull.enableAddedMassMoment = true;
    seedParamsFull.enableLEV             = false;
    seedParamsFull.levApplicationPoint   = 'colocated';

    % --- LEV constants for THIS seed (used only when enableLEV is true) -------
    % Every number below is derived or chosen as follows (full provenance in
    % levPlanformConstants and computeLEVForce):
    %   AR      = total span / mean chord -- "the ratio between the span and chord
    %             of the wing", as Rezgui et al. (2020) define it. The strip widths
    %             dz are arc lengths, so AR is unchanged by curvature. Test seed: 3.33.
    %   Kp, Ki  = Helmbold (1942) lift slope and Prandtl elliptic induced-drag
    %             factor for that AR -> 3.557 and 0.0955 for the test seed.
    %   Kv      = Kp - Kp^2*Ki, Rezgui et al. eq. (4) -> 2.349 for the test seed.
    %   RoCrit  = 3   anchor for the Rossby roll-off. NOT a published threshold:
    %                 Lentink & Dickinson (2009) show a stable LEV at Rg/c = 2.9 and
    %                 note real wings cluster near tip-radius Ro ~ 3. Tunable.
    %   p       = 4   sharpness of that roll-off. No literature value. Tunable.
    %   lambdaV = 0.5 fraction of the half-chord toward the leading edge, used only
    %                 by the 'forward' application point. No direct citation. Tunable.
    %   rossbyDefinition = 'kinematic' -- WHICH Rossby number the gate reads:
    %                 'kinematic' Ro = |v_ip|/(Omega c), Rossby's own ratio, with
    %                     Omega = |omega x s_hat| the strip's revolution rate. Reduces
    %                     exactly to r/c under pure revolution and shuts the gate for
    %                     a seed that is not revolving at all. The default.
    %                 'geometric' Ro = r/c, that special case applied as geometry.
    %                     Kept so pre-kinematic runs stay reproducible; it grants a
    %                     stable LEV to a parachuting or fluttering seed, which is
    %                     the defect the kinematic form fixes.
    % Override any of them with cfg.lev (e.g. cfg.lev.Kp, cfg.lev.RoCrit);
    % buildSeedParams keeps Kv consistent with eq. (4) when Kp or Ki change.
    stripChord = seedParamsFull.strips.chord;   stripDz = seedParamsFull.strips.dz;
    levSpan    = sum(stripDz);
    levChord   = sum(stripChord .* stripDz) / levSpan;
    lev        = levPlanformConstants(levSpan / levChord);
    lev.RoCrit  = 3;
    lev.p       = 4;
    lev.lambdaV = 0.5;
    lev.rossbyDefinition = 'kinematic';
    seedParamsFull.lev = lev;

    % --- Drop the switches the 3D RHS no longer honours --------------------
    % setupSeedShapeAndMass stamps the full planar switch set. seed6DOFODE3D no
    % longer reads these -- the terms they gated were removed:
    %   Sep-2026 audit : the span torque with its migrating CoP and reduced-
    %                    frequency attenuation; the Tx / Ty spin-damping torques.
    %   phase 4        : the span force itself (measured inert to <= 7e-4 relative
    %                    on every metric) and the velocity-sampling switch that
    %                    only it used. Edge drag replaces it.
    % Leaving them on the struct would advertise knobs that silently do nothing, so
    % strip them: a shape3d seedParams carrying one is a stale config, and now looks
    % like one. The planar model still has all of them.
    deadSwitches = {'enableSpanTorque', 'enableSpanCOPMigration', ...
                    'enableSpanTorqueAttenuation', 'enableTxDamping', ...
                    'enableNormalSpinDamping', 'enableSpanForce', ...
                    'enableSpanGeomVelocity'};
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
