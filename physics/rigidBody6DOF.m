function [xdot, core] = rigidBody6DOF(x, mp, F_aero_body, tau_body, gMag)
% RIGIDBODY6DOF  Shape-agnostic 6-DOF rigid-body state derivative.
%
% The part of the seed equations of motion that does NOT depend on the seed's
% shape. Given the current state, the mass properties (including added mass), and
% the NET aerodynamic force/torque already summed in the BODY frame, it adds
% gravity, solves translation and rotation with the added-mass terms, integrates
% the quaternion, and assembles the 13-state derivative.
%
% seed6DOFODE3D computes (F_aero_body, tau_body) and calls this core. (The archived
% planar seed6DOFODE fed it the same way; its switches are absent, so it ran the
% uncoupled path below unchanged.)
%
% INPUTS
%   x           : 13x1 state [ r(3); q(4); v(3); omega(3) ].
%                   q     body->world quaternion, scalar-first (normalised here)
%                   v     CoM velocity, INERTIAL frame (m/s)
%                   omega angular velocity, BODY frame (rad/s)
%                 (r is not used for the derivative.)
%   mp          : mass-properties struct from getMassProperties: .M .c .I_G
%                 .I_G_dot .A_trans .A_rot .A_coup and the switches below, plus
%                 .tau_inSectionMunk (set by seed6DOFODE3D; see MUNK MOMENT).
%   F_aero_body : 3x1 net aerodynamic force,  BODY frame (N).
%   tau_body    : 3x1 net aerodynamic torque, BODY frame (N*m), about the CoM.
%   gMag        : gravity magnitude (m/s^2); gravity acts along world -Y (Y-up).
%
% OUTPUTS
%   xdot : 13x1 state derivative [ v; dq; a_inertial; alpha_body ].
%   core : (optional) intermediates for post-processing:
%          .F_aero_inertial, .F_total_inertial (aero + gravity),
%          .F_addedMassRate (inertial: the added-mass velocity force, zero if off),
%          .tau_addedMass   (body: the added-mass velocity moment REMOVED from
%                            tau, i.e. tau_solve = tau - tau_addedMass),
%          .a_inertial, .alpha_body.
%
% =========================================================================
% THE ADDED-MASS PHYSICS: Kirchhoff's equations (Lamb, Hydrodynamics ch. VI)
% =========================================================================
% The fluid's kinetic energy, built from the strips by getAddedMass, is
%     T_A = 1/2 [V; W]' [A_trans  A_coup ; A_coup'  A_rot] [V; W]     (about the CoM)
% With M6 = blkdiag(M*I, I_G) + that matrix, P = (M + A_trans)V + A_coup W and
% L = A_coup' V + (I_G + A_rot) W, the body-frame equations are
%     M6 [Vdot; Wdot] = [ F - W x P ;  tau - W x L - V x P - I_G_dot W ]
% Every switch below selects pieces of this one system.
%
% ADDED-MASS RATE (mp.enableAddedMassRate; shape3d default ON)
% The W x (added part of P) force: spinning changes the entrained fluid's momentum
% even at constant speed. Without coupling it is solved in the INERTIAL frame as
%     (M*I + A_i)*a = F - Adot_i*v,   Adot_i*v = w x (A_i*v) - A_i*(w x v)
% (A_i = R*A_trans*R'; verified against a finite difference of R*A_b*R' to 6.2e-9,
% and measured at 31-43% of the net force in spinning cases). It is the
% inertial-frame form of minimal_imp's (m+m2)*omega*vyp and -(m+m1)*omega*vxp
% (testNewPhysicsTerms A5).
%
% MUNK MOMENT (mp.enableAddedMassMoment, shape3d default ON;
%              mp.addedMassMomentMode = 'residual' (shape3d default) | 'full')
% The V x P term: the moment that turns a body broadside-on to its own motion.
%   'full'     : all of it, -V x (A V). For a flat plate (m2 v w, 0, -m2 u v).
%   'residual' : minus the IN-SECTION part, the 2D Munk moment of each strip about
%                its own span axis at its own centre velocity,
%                    tau_IS = -sum_i m2_i vc_i vn_i s_i      (seed6DOFODE3D)
%                That part is already inside the strips' measured centre of
%                pressure: minimal_imp's l_cp is the TOTAL centre of pressure (0.298
%                chord ahead of centre at alpha -> 0, where potential flow puts the
%                whole force at the quarter chord), and its small-angle moment slope
%                -1.564/rad is the Munk moment's -pi/2. Adding the Munk pitch moment
%                on top doubled the pitching moment at small angles and gave 6-7x it
%                past stall. What remains is the SPANWISE-flow part (roll, from
%                normal velocity with spanwise sliding), which the strips cannot see
%                and which stops a plate falling span-first. In 2D 'residual'
%                removes the whole term, as minimal_imp does.
%   (The planar-era note here quoted APW eq. (6.3) with both a Munk term and a
%    centre-of-pressure torque; that is not compatible with minimal_imp's l_cp.)
%
% OFF-CENTRE COUPLING (mp.enableAddedMassCoupling, shape3d default ON;
%                      mp.addedMassCouplingForm = 'symmetric' (default) | 'oneSided')
% The added mass is attached at the strip centres, not at the CoM, so translation
% and rotation couple through A_coup. ON solves the 6x6 system above in the body
% frame. In 2D (CoM lcm ahead of the plate centre) it gives
%     x':  -m2 omega^2 lcm          y':  +m2 lcm omegadot      (minimal_imp has both)
%     torque: +m2 lcm vdot                                      (minimal_imp lacks it)
% 'symmetric' is Kirchhoff's matrix: energy is conserved in ideal flow (5e-12 over
% 2 s). 'oneSided' drops the A_coup' block of the ROTATION rows, reproducing
% minimal_imp's equations exactly; it exists for testMinimalImpEquivalence and is
% not physics (it drifts 12% in energy over 2 s in ideal flow).
% OFF keeps the uncoupled inertial-frame solve, byte-identical to before.

    % --- Parse state (normalise q defensively; state may have drifted) ------
    q     = x(4:7);   q = q / norm(q);
    v     = x(8:10);          % CoM velocity, inertial
    omega = x(11:13);         % angular velocity, body frame
    R     = quatToRotm(q);    % body->world

    % --- Forces: aero (body->inertial) + gravity ---------------------------
    F_aero_inertial  = R * F_aero_body;
    F_grav           = mp.M * [0; -gMag; 0];        % world -Y; buoyancy neglected
    F_total_inertial = F_aero_inertial + F_grav;

    useRate   = isfield(mp, 'enableAddedMassRate')   && mp.enableAddedMassRate;
    useMoment = isfield(mp, 'enableAddedMassMoment') && mp.enableAddedMassMoment;
    residual  = useMoment && isfield(mp, 'addedMassMomentMode') && ...
                strcmp(mp.addedMassMomentMode, 'residual');
    useCoup   = isfield(mp, 'enableAddedMassCoupling') && mp.enableAddedMassCoupling;

    if useCoup
        % ===== Coupled path: Kirchhoff's 6x6 system in the body frame ========
        V   = R.' * v;
        Fb  = R.' * F_total_inertial;
        At  = mp.A_trans;   Ac = mp.A_coup;   Ar = mp.A_rot;
        PA  = At*V + Ac*omega;                       % fluid momentum
        fF  = Fb - cross(omega, mp.M*V);             % body rotating-frame term
        fT  = tau_body - cross(omega, (mp.I_G + Ar)*omega) - mp.I_G_dot*omega;
        F_rate_b = [0; 0; 0];   tau_addedMass = [0; 0; 0];
        if useRate
            F_rate_b = cross(omega, PA);
            fF = fF - F_rate_b;
            fT = fT - cross(omega, Ac.'*V);
        end
        if useMoment
            tau_addedMass = cross(V, PA);
            if residual; tau_addedMass = tau_addedMass + mp.tau_inSectionMunk; end
            fT = fT - tau_addedMass;
        end
        if isfield(mp, 'addedMassCouplingForm') && strcmp(mp.addedMassCouplingForm, 'oneSided')
            lowerLeft = zeros(3);
        else
            lowerLeft = Ac.';
        end
        M6  = [mp.M*eye(3) + At, Ac;  lowerLeft, mp.I_G + Ar];
        acc = M6 \ [fF; fT];
        a_inertial = R * (acc(1:3) + cross(omega, V));
        alpha_body = acc(4:6);
        F_addedMassRate = R * F_rate_b;
    else
        % ===== Uncoupled path (no A_coup): unchanged since phase 5 ============
        % --- Added-mass rate:  F_solve = F_total - Adot_i*v -----------------
        F_addedMassRate = [0; 0; 0];
        if useRate
            A_i     = R * mp.A_trans * R.';              % added mass, inertial frame
            omega_i = R * omega;                         % angular velocity, inertial
            F_addedMassRate = cross(omega_i, A_i*v) - A_i*cross(omega_i, v);   % = Adot_i*v
        end
        F_solve = F_total_inertial - F_addedMassRate;

        % --- Munk moment:  tau_solve = tau - v x (A*v) [- tau_IS] -----------
        % Formed in the BODY frame, where A_trans lives and where the torque is
        % applied, so no rotation of the tensor is needed.
        tau_addedMass = [0; 0; 0];
        if useMoment
            v_body        = R.' * v;
            tau_addedMass = cross(v_body, mp.A_trans * v_body);      % = v x (A v)
            if residual; tau_addedMass = tau_addedMass + mp.tau_inSectionMunk; end
        end
        tau_solve = tau_body - tau_addedMass;

        a_inertial = translationDynamics(F_solve, R, mp);            % with added mass
        alpha_body = rotationDynamics(tau_solve, omega, mp);         % modified Euler
    end

    % --- Orientation kinematics --------------------------------------------
    dq = quatKinematics(q, omega);

    % --- Assemble (dr/dt = v, since inertial velocity is a state) -----------
    xdot = [ v ; dq ; a_inertial ; alpha_body ];

    if nargout > 1
        core.F_aero_inertial  = F_aero_inertial;
        core.F_total_inertial = F_total_inertial;
        core.F_addedMassRate  = F_addedMassRate;
        core.tau_addedMass    = tau_addedMass;
        core.a_inertial       = a_inertial;
        core.alpha_body       = alpha_body;
    end
end
