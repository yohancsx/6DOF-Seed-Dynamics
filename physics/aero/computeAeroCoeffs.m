function coeffs = computeAeroCoeffs(alpha, aero)
% COMPUTEAEROCOEFFS  Quasi-steady aerodynamic coefficients vs angle of attack.
%
% Direct extraction of the lift/drag/center-of-pressure and rotational
% coefficient mathematics from the 2D model (minimal_imp_Commented.m, the
% Andersen-Pesavento-Wang style flat-plate model). This returns COEFFICIENTS
% ONLY -- the dimensional force/torque assembly (multiplying by 0.5*rho*chord,
% speed, omega, moment arms, etc.) is done elsewhere.
%
% The 2D code computed everything per unit depth for a single plate of chord
% l; here the chord-dependence is removed so the coefficients are reusable per
% strip. The center of pressure is returned as a FRACTION of chord; multiply by
% the local strip chord downstream to get a physical length.
%
% INPUTS
%   alpha : local angle of attack (rad). Scalar OR array (vectorised over
%           strips). DEFINITION must match the 2D model: alpha is the angle of
%           the strip's in-plane RELATIVE VELOCITY measured from the body
%           x'-axis (chord), i.e. alpha = atan2(v_normal, v_chord). In the 2D
%           code this was atan2(vyp - omega*lcm, vxp). Compute alpha from the
%           strip velocity (not the negated wind) so the branch logic lines up.
%   aero  : (optional) struct of empirical constants. If omitted, the exact
%           minimal_imp values are used (see defaultAeroParams below).
%           Fields: CL1, CL2, CD0, CD1, CD2, CR, CCP0, CCP1, CCP2, as, d, C_fy,
%           C_span.
%
% OUTPUT (struct, each field matches the size of alpha unless noted)
%   coeffs.CT        : lift / normal-force coefficient (the 2D "CT_alpha").
%                      alpha-dependent.
%   coeffs.CD        : drag coefficient (the 2D "CD_alpha"). alpha-dependent.
%   coeffs.l_cp_frac : center of pressure as a FRACTION of chord, measured from
%                      the geometric center (the 2D "l_cp" before the *l).
%                      alpha-dependent.
%   coeffs.CR        : rotational-lift (added-circulation) coefficient. The
%                      coefficient on the "-l*CR*omega" term in L_0 and on the
%                      Tcr rotational-lift torque. CONSTANT (alpha-independent).
%   coeffs.CD_rot    : rotational spin-damping coefficient (= CD2). The
%                      coefficient on the Tr ~ omega*|omega| damping torque.
%                      CONSTANT (alpha-independent).
%   coeffs.CD0       : zero-incidence/parasitic drag coefficient. Used
%                      directly as the spin-damping coefficient about the
%                      plate-normal (body y) axis (normalSpinDamping.m).
%                      CONSTANT (alpha-independent).
%   coeffs.C_fy      : tuning factor for the plate-normal spin-damping torque
%                      (normalSpinDamping.m); no minimal_imp analogue.
%                      CONSTANT (alpha-independent).
%   coeffs.C_span    : tuning factor for the whole-seed spanwise-flow force
%                      (computeSpanForce.m); no minimal_imp analogue.
%                      CONSTANT (alpha-independent).

% -------------------------------------------------------------------------
% 0. PARAMETERS
% -------------------------------------------------------------------------
if nargin < 2 || isempty(aero)
    aero = defaultAeroParams();
end

CL1 = aero.CL1;  CL2 = aero.CL2;                       % lift constants
CD0 = aero.CD0;  CD1 = aero.CD1;  CD2 = aero.CD2;      % drag constants
CR  = aero.CR;                                         % rotational-lift constant
CCP0 = aero.CCP0; CCP1 = aero.CCP1; CCP2 = aero.CCP2;  % CoP constants
as = aero.as;    d = aero.d;                           % stall angle, blend width
C_fy = aero.C_fy;                                      % y-axis spin-damping tuning factor
C_span = aero.C_span;                                  % spanwise-flow force tuning factor

% -------------------------------------------------------------------------
% 1. SMOOTH ATTACHED<->SEPARATED BLEND WEIGHTS  (2D lines 122-123)
%    tanh blends the attached model (small |alpha|) into the separated model
%    (large |alpha|) around stall angle "as" with width "d".
% -------------------------------------------------------------------------
absA = abs(alpha);
f1 = ( 1 - tanh((absA - as)/d) )/2;        % ~1 attached, ~0 stalled (|alpha|<pi/2)
f2 = ( 1 - tanh((pi - absA - as)/d) )/2;   % same blend mirrored near alpha=+/-pi

% -------------------------------------------------------------------------
% 2. LIFT / NORMAL-FORCE COEFFICIENT  (2D lines 129-131)
%    Blend of attached (CL1*sin alpha) and separated (CL2*sin 2alpha) lift.
%    Three branches; sign flips set the correct lift direction front/back.
% -------------------------------------------------------------------------
CT_1 =  CL1.*f1.*sin(alpha)      + (1-f1).*CL2.*sin(2*alpha);                  % |alpha| < 90 deg
CT_2 = -( CL1.*f2.*(pi-absA)     + (1-f2).*CL2.*sin(2*(pi-absA)) );            % alpha in ( 90,180)
CT_3 =    CL1.*f2.*(pi-absA)     + (1-f2).*CL2.*sin(2*(pi-absA));              % alpha in (-180,-90)
% NOTE (preserved from 2D): CT_2/CT_3 use (pi-abs(alpha)) directly where CT_1
% uses sin(alpha) -- a small-angle approximation of sin(pi-|alpha|) near
% alpha=+/-pi. Kept identical to the 2D code; flag if you intend to change it.

% -------------------------------------------------------------------------
% 3. DRAG COEFFICIENT  (2D lines 139-140)
%    Attached: CD0 + CD1*sin^2 alpha ; separated: CD2*sin^2 alpha. Blended.
% -------------------------------------------------------------------------
CD_1 = f1.*(CD0 + CD1.*sin(alpha).^2)    + (1-f1).*CD2.*sin(alpha).^2;         % |alpha| < 90
CD_2 = f2.*(CD0 + CD1.*(pi-absA).^2)     + (1-f2).*CD2.*sin(pi-alpha).^2;      % |alpha| > 90
% NOTE (preserved from 2D): CD_2 uses (pi-abs(alpha))^2 in the first term (vs
% sin^2 in CD_1) and sin(pi-alpha) (no abs) in the second. The squaring makes
% the missing abs harmless, but the angle-vs-sine choice mirrors the CT note.

% -------------------------------------------------------------------------
% 4. CENTER OF PRESSURE (fraction of chord from geometric center)  (2D 146-152)
%    Attached branch: quadratic-in-angle; separated branch: linear ramp.
% -------------------------------------------------------------------------
r  = absA;        % |alpha|
rr = pi - r;      % supplement of |alpha| (high-angle branch)
lcp_1 =   f1.*(CCP0 - CCP1.*r.^2)  + (1-f1).*(CCP2 - (CCP2/(pi/2)).*r);        % low-angle fraction
lcp_2 = -( f2.*(CCP0 - CCP1.*rr.^2) + (1-f2).*(CCP2 - (CCP2/(pi/2)).*rr) );    % high-angle fraction (sign-flipped)

% -------------------------------------------------------------------------
% 5. BRANCH SELECTION  (exact replica of the 2D if/elseif/else, lines 155-167)
%    Vectorised with logical masks so alpha may be scalar or an array:
%      m1 : |alpha| <= 90 deg          -> attached/forward  (CT_1, CD_1, lcp_1)
%      m2 : alpha > 90 deg             -> reversed, front    (CT_2, CD_2, lcp_2)
%      m3 : alpha < -90 deg            -> reversed, back     (CT_3, CD_2, lcp_2)
%    Drag and CoP are identical for m2 and m3 (only the lift sign differs).
% -------------------------------------------------------------------------
m1 = (r <= pi/2);        % attached/forward
m2 = (alpha > pi/2);     % reversed front
m3 = ~(m1 | m2);         % reversed back (alpha < -pi/2)

CT  = CT_1.*m1  + CT_2.*m2 + CT_3.*m3;
CD  = CD_1.*m1  + CD_2.*(m2 | m3);
lcp = lcp_1.*m1 + lcp_2.*(m2 | m3);

% -------------------------------------------------------------------------
% 6. PACK OUTPUT
% -------------------------------------------------------------------------
coeffs.CT        = CT;     % lift/normal-force coefficient (alpha-dependent)
coeffs.CD        = CD;     % drag coefficient              (alpha-dependent)
coeffs.l_cp_frac = lcp;    % CoP fraction of chord         (alpha-dependent)
coeffs.CR        = CR;     % rotational-lift coefficient   (constant)
coeffs.CD_rot    = CD2;    % rotational damping coefficient (constant, = CD2)
coeffs.CD0       = CD0;    % zero-incidence drag coefficient (constant)
coeffs.C_fy      = C_fy;   % y-axis spin-damping tuning factor (constant)
coeffs.C_span    = C_span; % spanwise-flow force tuning factor (constant)

end   % computeAeroCoeffs


% =========================================================================
% LOCAL: exact empirical constants from minimal_imp_Commented.m (lines 81-89,
%        120-121). CGS-model fit values; dimensionless except as,d (radians).
% =========================================================================
function aero = defaultAeroParams()
    aero.CL1  = 5.2;          % attached-flow lift slope
    aero.CL2  = 0.95;         % separated-flow lift amplitude (sin 2alpha term)
    aero.CD0  = 0.1;          % parasitic (zero-incidence) drag
    aero.CD1  = 5;            % attached-flow induced drag (~sin^2 alpha)
    aero.CD2  = 1.9;          % separated/broadside drag (also spin damping)
    aero.CR   = 1.1;          % rotational-lift / added-circulation coefficient
    aero.CCP0 = 0.3;          % CoP: baseline offset (attached)
    aero.CCP1 = 3.5;          % CoP: curvature vs angle (attached)
    aero.CCP2 = 0.2;          % CoP: high-angle offset (separated)
    aero.as   = 14*pi/180;    % stall angle (~14 deg)
    aero.d    =  6*pi/180;    % blend width (~6 deg)
    aero.C_fy = 1.0;          % y-axis (plate-normal) spin-damping tuning factor
                               % (no minimal_imp analogue; see normalSpinDamping.m)
    aero.C_span = 1.0;        % spanwise-flow force tuning factor
                               % (no minimal_imp analogue; see computeSpanForce.m)
end