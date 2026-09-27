function [F_lev_local, xApp, info] = computeLEVForce(vChord, vNormal, alpha, chord, dz, r, Omega, lcpFrac, lev, appPoint, rhoFluid)
% COMPUTELEVFORCE  Leading-edge-vortex (LEV) vortex lift on one strip, strip frame.
%
% =========================================================================
% PHYSICS, PLAINLY
% =========================================================================
% A plate at incidence "wants" a suction peak right at its leading edge. With a
% sharp edge the flow separates there instead, and that suction is lost. If the
% separated shear layer rolls up into a vortex that STAYS over the plate, the lost
% suction reappears as an extra force pressing on the surface (Polhamus' leading-
% edge suction analogy). On a revolving wing the vortex does stay put, because the
% rotation's centripetal and Coriolis accelerations hold it on -- which is how
% samaras get their unexpectedly high lift (Lentink et al. 2009, Science
% 324:1438-1440). On a plate that is merely translating, the vortex is shed.
%
% =========================================================================
% THE EQUATIONS AND WHERE EACH COMES FROM
% =========================================================================
%   (1) Vortex-lift magnitude
%         C_L,v = K_v * sin^2(a_e) * cos(a_e)
%       Rezgui, Arroyo & Theunissen (2020), Aeronautical J. 124(1278):1236-1261,
%       eq. (3) with sweep Lambda = 0; from Polhamus (1966, NASA TN D-3767; 1971,
%       J. Aircraft 8(4):193-199). a_e is the angle between the flow and the plate,
%       0..90 deg. Polhamus and Rezgui write the law for 0 <= alpha <= 90 deg only.
%       This model's alpha spans -180..180, so a_e is the SAME effective angle the
%       APW branches in computeAeroCoeffs use: a_e = |alpha| up to 90 deg and
%       pi - |alpha| beyond. Written in alpha directly that is
%       sin^2(alpha)*|cos(alpha)|, which is what is coded.
%       ONLY this vortex term is added. Rezgui et al.'s potential term,
%       K_p sin(a) cos^2(a) (their eq. 1), is NOT -- the existing APW coefficient
%       law already supplies the potential/separated lift, and adding eq. (1) too
%       would count it twice.
%
%   (2) Sign -- taken from the APW branch, never hand-derived
%         |alpha| <= 90 deg : sign(alpha)   (as CT_1 ~ sin(alpha), sin(2 alpha))
%         alpha   >  90 deg : -1            (as CT_2)
%         alpha   < -90 deg : +1            (as CT_3)
%       so the vortex lift always pushes the same way the APW lift does.
%
%   (3) Rossby gate -- how much of the vortex lift this strip is allowed
%         G = 1 / (1 + (Ro/Ro_crit)^p)
%       Basis: Lentink & Dickinson (2009), J. Exp. Biol. 212:2705-2719 -- the LEV
%       is held on by the rotation's centripetal and Coriolis accelerations, which
%       scale as 1/Ro. They observe a stable LEV on a revolving wing at
%       Ro = Rg/c = 2.9, an unstable one on a translating wing (Ro = inf), and force
%       coefficients that change with Ro over 2.9 -> 3.6 -> 4.4. They do NOT give
%       a critical Rossby number or a transfer function. THIS GATE IS OUR
%       CONSTRUCTION: Ro_crit (default 3) anchors the roll-off at that scale, and
%       p (default 4) sets how sharp it is. Both are tuning parameters. Rezgui et
%       al. apply the vortex lift everywhere, ungated; G = 1 reproduces that.
%
%       WHICH Ro -- lev.rossbyDefinition, a switch:
%
%       'kinematic' (DEFAULT)      Ro = |v_ip| / (Omega * c)
%          Rossby's own ratio (Rossby 1936): the speed the section actually sees,
%          over the rotation rate times the length scale. It is the quantity
%          Lentink & Dickinson are forming when they scale the Coriolis and
%          centripetal accelerations "with respect to the fluid's convective
%          acceleration"; their Ro = Rg/c is its PURE-REVOLUTION special case.
%          Omega is the rate at which this strip REVOLVES about the CoM,
%          |omega x s_hat| -- rotation about the strip's own span axis is pitching,
%          not revolution, so it is excluded. Under pure revolution v_ip = Omega*r,
%          so this reduces EXACTLY to r/c (testNewPhysicsTerms L3b, to ~1e-16),
%          and a merely translating plate (Omega -> 0) gets Ro -> inf, G -> 0: no
%          stable vortex, which is what is observed.
%          Coded DIVISION-FREE, algebraically the same function:
%              G = (Ro_crit*c*Omega)^p / ( (Ro_crit*c*Omega)^p + |v_ip|^p )
%          so nothing divides by Omega and the gate stays finite and smooth for the
%          integrator as the seed stops revolving. With neither rotation nor flow
%          (0/0) G is set to 0; the force is zero there regardless.
%
%       'geometric'                Ro = r / c
%          The pure-revolution special case used as a geometric property of the
%          strip. Kept so runs made before the kinematic gate existed remain
%          reproducible. ITS KNOWN FLAW is why the default moved: r > 0 whether or
%          not the seed revolves, so it grants a stable LEV to a parachuting,
%          gliding or fluttering seed, which cannot hold one.
%
%   (4) Force -- the same lift convention as computeStripForces
%         F = 1/2 * rho * c * dz * (G * C_L,v) * |v| * [v_n ; -v_c ; 0]
%       i.e. perpendicular to the strip's in-plane velocity, exactly as the APW
%       translational lift is formed. It is lift, so it does no work in the
%       strip plane. It is NOT scaled by the per-strip liftMult.
%
%   (5) Application point -- a SWITCH, because the literature conflicts
%       'colocated' (default): at the existing APW centre of pressure,
%            x_app = l_cp_frac * c
%          Snyder & Lamar (1972), NASA TN D-6994: on low-aspect-ratio delta wings
%          the vortex-lift load distribution is "similar in shape to that of the
%          potential-flow longitudinal loading", i.e. it acts at about the same
%          centre of pressure. With this choice the LEV adds lift but does NOT move
%          the centre of pressure.
%       'forward': toward the leading edge,
%            x_app = lambda_v * (c/2) * cos(alpha)
%          lambda_v is the fraction of the half-chord from mid-chord toward the
%          leading edge (0 = mid-chord, 1 = on the edge; default 0.5, i.e. about a
%          quarter-chord aft of the leading edge). NO DIRECT CITATION: motivated by
%          the LEV sitting over the front of an UNSWEPT section, which a delta wing
%          is not. cos(alpha) points it at the upwind (leading) edge -- sign(cos a)
%          = sign(v_chord), the same quantity the APW branches and l_cp track -- and
%          sends it to zero at broadside, where there is no distinct leading edge.
%          Since the force also vanishes there, the moment stays continuous as the
%          leading edge flips.
%       x_app is measured along the strip's chord axis from its geometric centre.
%
% =========================================================================
% DEFINITIONS FOR THIS MODEL
% =========================================================================
%   v_ip  The strip's IN-PLANE speed, hypot(v_chord, v_normal) -- the same velocity
%         the APW coefficients are evaluated at, so the gate responds to the flow
%         the section actually sees, descent included. (Descent along the spin axis
%         raises the convective acceleration without raising the Coriolis one, so it
%         belongs in the numerator: more through-flow, relatively weaker rotation.)
%   Omega The rate at which this strip REVOLVES about the CoM, |omega x s_hat|,
%         with omega the body angular velocity. Per strip, so a twisted or curved
%         seed is handled correctly. Excludes omega . s_hat, which pitches the
%         strip about its own span axis rather than carrying it around. The
%         narrower alternative |omega . n_hat| (only rotation that sweeps the strip
%         chordwise) is defensible too; |omega x s_hat| is used because it is the
%         full non-pitching rotation.
%   r     Distance of the strip from the seed's rotation axis, taken as its
%         SPANWISE distance from the CoM along the local span axis,
%         r = |s_hat . (p_strip - c_CoM)|  (flat seed: |z_strip - z_CoM|). Used only
%         by the 'geometric' definition. Computed by the caller every step, so it
%         tracks a moving CoM. Lentink & Dickinson use a single whole-wing radius
%         (radius of gyration, or tip radius in their Eq. 6); ours is per-strip, so
%         the gate varies along the span.
%   c     The strip's local chord. Lentink & Dickinson use the mean chord.
%   Ro    LOCAL Rossby number (per strip, per instant) -- not the same quantity as
%         their whole-wing Ro, though it reduces to it under pure revolution.
%
% INPUTS
%   vChord, vNormal : strip velocity components in its local (chord, normal) plane (m/s).
%   alpha           : angle of attack, atan2(vNormal, vChord) (rad), -pi..pi.
%   chord, dz       : strip chord and spanwise width (m).
%   r               : strip distance from the rotation axis (m)   -- 'geometric' only.
%   Omega           : strip revolution rate |omega x s_hat| (rad/s) -- 'kinematic' only.
%   lcpFrac         : APW centre-of-pressure fraction at this alpha (computeAeroCoeffs).
%   lev             : struct with .Kv, .RoCrit, .p, .lambdaV and .rossbyDefinition
%                     ('kinematic' default | 'geometric') -- see setupSeedShape3D.
%   appPoint        : 'colocated' | 'forward'.
%   rhoFluid        : fluid density (kg/m^3).
%
% OUTPUTS
%   F_lev_local : 3x1 vortex-lift force in the strip frame [chord; normal; span] (N).
%   xApp        : chordwise application offset from the strip centre (m).
%   info        : struct with .CT_lev (gated coefficient, signed), .G, .Ro, .Omega.
%                 .Ro is Inf where the strip is not revolving (consistent with G = 0).

    % --- (3) Rossby gate ----------------------------------------------------
    v_ip = hypot(vChord, vNormal);          % in-plane speed the section sees (m/s)

    if isfield(lev, 'rossbyDefinition')
        roDef = char(lev.rossbyDefinition);
    else
        roDef = 'kinematic';                % default for a struct built before the switch
    end

    switch roDef
        case 'kinematic'
            % G = (Ro_crit c Omega)^p / ((Ro_crit c Omega)^p + v_ip^p)
            %   == 1/(1 + (Ro/Ro_crit)^p) with Ro = v_ip/(Omega c), but with no
            %   division by Omega, so Omega -> 0 is finite and smooth.
            scale = lev.RoCrit * chord * Omega;
            num   = scale ^ lev.p;
            den   = num + v_ip ^ lev.p;
            if den > 0
                G = num / den;
            else
                G = 0;              % neither revolving nor moving: no vortex (F = 0 anyway)
            end
            if Omega > 0
                Ro = v_ip / (Omega * chord);
            else
                Ro = Inf;           % reported value consistent with G = 0
            end
        case 'geometric'
            Ro = r / chord;
            G  = 1 / (1 + (Ro / lev.RoCrit)^lev.p);
        otherwise
            error('computeLEVForce:badRossbyDefinition', ...
                  'lev.rossbyDefinition must be ''kinematic'' or ''geometric'' (got ''%s'').', roDef);
    end

    % --- (1)+(2) signed, gated vortex-lift coefficient ---------------------
    mag = lev.Kv * sin(alpha)^2 * abs(cos(alpha));      % K_v sin^2(a_e) cos(a_e)
    if abs(alpha) <= pi/2
        sgn = sign(alpha);                              % APW branch CT_1
    elseif alpha > pi/2
        sgn = -1;                                       % APW branch CT_2
    else
        sgn = 1;                                        % APW branch CT_3
    end
    CT_lev = G * sgn * mag;

    % --- (4) force, lift convention of computeStripForces ------------------
    Lt0  = 0.5 * rhoFluid * chord * dz * CT_lev * v_ip;
    F_lev_local = [ Lt0*vNormal ; -Lt0*vChord ; 0 ];

    % --- (5) application point ---------------------------------------------
    switch appPoint
        case 'colocated'
            xApp = lcpFrac * chord;
        case 'forward'
            xApp = lev.lambdaV * (chord/2) * cos(alpha);
        otherwise
            error('computeLEVForce:badAppPoint', ...
                  'levApplicationPoint must be ''colocated'' or ''forward'' (got ''%s'').', appPoint);
    end

    info = struct('CT_lev', CT_lev, 'G', G, 'Ro', Ro, 'Omega', Omega);
end
