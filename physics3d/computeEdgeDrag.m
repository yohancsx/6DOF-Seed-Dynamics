function [F_edge, tau_edge, vSpan] = computeEdgeDrag(vApply_body, rArm_body, spanAxis_body, A_edge, C_d_edge, rhoFluid)
% COMPUTEEDGEDRAG  Crossflow drag on a seed sliding along its own span.
%
% PHYSICS, PLAINLY
% The strip model resolves flow in each strip's chord-normal plane only, so a
% seed sliding sideways along its span feels EXACTLY ZERO aerodynamic force --
% the strips produce none by construction. Physically it presents its wingtip to
% the flow: a small bluff cross-section, thickness by chord. This is the smallest
% honest term that removes that unphysical zero:
%
%     F_edge = -1/2 * rho * C_d * A_edge * |v_s| * v_s * s_hat
%     A_edge = thickness * mean chord
%
% where v_s is the velocity component along the span axis s_hat. It is PURE drag:
% it opposes the spanwise velocity and nothing else. No lift, no migrating centre
% of pressure, no separate torque coefficient.
%
% WHY THE APPLICATION POINT IS THE AREA CENTROID
% The force is directed along the span, so WHERE along the span it acts does not
% change its moment -- sliding the point along the force's own line of action
% leaves r x F unchanged. Only the force line's chordwise/normal position matters,
% and for a chord-symmetric section that line runs through the area centroid.
% So applying at the centroid is exact for the moment, and the moment is just
% cross(r_arm, F).
%
% WHAT IT REPLACES
% computeSpanForce, which applied the full 2D lift+drag law in the span-normal
% plane with the span in the chord's role, scaled by a tuned C_span = 0.2. Phase 3
% measured it inert (<= 7e-4 relative on every benchmark metric), and it had the
% wrong shape: smallest for a pure spanwise slide (where it borrowed the attached
% CD0) and largest where the strips already supply the force. Expect THIS term to
% be inert too -- it is ~4x smaller again. It is here to remove a known-unphysical
% zero, not to move results.
%
% KNOWN OMISSION
% Only the tip form drag is modelled. Laminar skin friction along the two faces is
% of comparable size at this Reynolds number (Re_span ~ 3e3, Cf ~ 1.33/sqrt(Re)
% gives ~2e-5 N at 1 m/s, vs ~2.2e-5 N form drag), so this is a LOWER BOUND on the
% true edgewise resistance, by roughly a factor of two.
%
% INPUTS
%   vApply_body   : 3x1 velocity of the application point (the area centroid),
%                   BODY frame, INCLUDING the omega x r transport (m/s).
%   rArm_body     : 3x1 application point relative to the CoM, BODY frame (m).
%   spanAxis_body : 3x1 unit span axis, BODY frame ([0;0;1] for this layout).
%   A_edge        : frontal area presented to spanwise flow, thickness*chord (m^2).
%   C_d_edge      : bluff-tip drag coefficient (~1.2 for a rectangular edge).
%   rhoFluid      : fluid density (kg/m^3).
%
% OUTPUTS
%   F_edge   : 3x1 edge-drag force, BODY frame (N).
%   tau_edge : 3x1 its moment about the CoM, BODY frame (N*m).
%   vSpan    : scalar spanwise velocity component used (m/s), for diagnostics.

    vSpan    = spanAxis_body.' * vApply_body;                          % v_s
    F_edge   = -0.5 * rhoFluid * C_d_edge * A_edge * abs(vSpan) * vSpan ...
               * spanAxis_body;                                        % opposes v_s
    tau_edge = cross(rArm_body, F_edge);
end
