function lev = levPlanformConstants(AR)
% LEVPLANFORMCONSTANTS  Polhamus vortex-lift constants for a flat wing of aspect ratio AR.
%
% The leading-edge-vortex (LEV) term in this model follows Rezgui, Arroyo &
% Theunissen (2020), "Model for sectional leading-edge vortex lift for the
% prediction of rotating samara seeds performance", The Aeronautical Journal
% 124(1278):1236-1261, doi:10.1017/aer.2020.25. They adapt Polhamus' leading-edge
% suction analogy (Polhamus 1966, NASA TN D-3767; Polhamus 1971, J. Aircraft
% 8(4):193-199, doi:10.2514/3.44254) to a samara blade section. Their vortex-lift
% coefficient needs ONE constant, K_v, and it is NOT a free tuning parameter --
% it follows from the planform:
%
%   K_v = K_p - K_p^2 * K_i                       Rezgui et al. eq. (4)
%
%   K_p : the wing's lift-curve slope (per radian), dC_L/d(alpha) at small alpha
%         (their eq. 2).
%   K_i : the induced-drag factor, K_i = d(C_Di) / d(C_L^2)  (their eq. 5).
%
% Rezgui et al. obtain K_p and K_i from a vortex-lattice code (Tornado) for their
% planform. No vortex-lattice code is part of this project, so both are estimated
% here from standard finite-wing theory for a flat, unswept wing. The two MUST come
% from the same wing model: mixing this model's 2D APW slope (CL1 = 5.2) with a 3D
% K_i is inconsistent and, at low aspect ratio, spuriously drives K_v toward zero.
%
% HOW EACH NUMBER IS CALCULATED
%   a0  = 2*pi         Thin-aerofoil lift slope of an infinite flat plate, per rad.
%   K_p = a0 / ( sqrt(1 + (a0/(pi*AR))^2) + a0/(pi*AR) )
%                      Helmbold (1942), "Der unverwundene Ellipsenflugel als
%                      tragende Wirbelflache" -- the standard lift slope of a
%                      straight wing of LOW aspect ratio. It tends to a0 as AR -> inf
%                      and to the slender-wing limit pi*AR/2 as AR -> 0, which is why
%                      it is preferred over Prandtl's a0/(1 + a0/(pi*AR)) here:
%                      Rezgui et al. note Prandtl's lifting line gives good results
%                      only above AR ~ 3, and this seed sits at or below that.
%   K_i = 1/(pi*AR)    Prandtl lifting-line induced drag with elliptic loading,
%                      C_Di = C_L^2/(pi*AR*e), span efficiency e = 1.
%   K_v = K_p - K_p^2*K_i
%                      Rezgui et al. (2020) eq. (4). With the two estimates above it
%                      is always positive, since Helmbold's K_p < pi*AR for every AR.
%
% FOR THIS PROJECT'S SEED (span 0.050 m, chord 0.015 m)
%   AR is the ratio of span to chord of the wing, as Rezgui et al. define it; for
%   a plate of varying chord the builder passes total span / mean chord. That gives
%   AR = 3.33 -> K_p = 3.557, K_i = 0.0955, K_v = 2.349.
%   Other readings of "the wing" give other values (testing/shape3d/exploreLEVTerm
%   prints them): a single half-span blade measured from a centred CoM, AR = 1.67,
%   gives K_v = 1.287; Rezgui et al.'s real samara, AR = 4.38, gives K_v = 2.853.
%
% INPUT
%   AR  : aspect ratio (span / mean chord), > 0.
% OUTPUT (struct)
%   lev.AR, lev.a0, lev.Kp, lev.Ki, lev.Kv   as defined above.

    if ~(isscalar(AR) && isfinite(AR) && AR > 0)
        error('levPlanformConstants:badAR', 'AR must be a positive finite scalar (got %g).', AR);
    end
    a0 = 2*pi;                                   % thin-aerofoil slope, per rad
    x  = a0 / (pi*AR);                           % = 2/AR
    Kp = a0 / (sqrt(1 + x^2) + x);               % Helmbold (1942)
    Ki = 1 / (pi*AR);                            % Prandtl lifting line, e = 1
    Kv = Kp - Kp^2*Ki;                           % Rezgui et al. (2020) eq. (4)

    lev = struct('AR',AR, 'a0',a0, 'Kp',Kp, 'Ki',Ki, 'Kv',Kv);
end
