function pc = planarAeroConstants(aero)
% PLANARAEROCONSTANTS  The archived planar model's tuning constants.
%
% These five constants scaled the planar model's empirical whole-seed terms (the
% span force and torque, and the Tx / Ty spin-damping torques). They lived in the
% shared computeAeroCoeffs until the planar model was archived (2026-10-15); the
% active model never read them, so they moved here with the code that does.
% Values are the planar builder-default working tuning, unchanged.
%
% The original, untouched planar model is at git tag planar-final.
%
% INPUT
%   aero : (optional) seedParams.aero override struct; any of the five fields
%          given there wins, exactly as the old defaults-merge did.
% OUTPUT
%   pc   : struct with .C_fy, .C_span, .C_span_torque, .k0_spanTorque, .C_Tx.

    pc.C_fy          = 1;     % y-axis (plate-normal) spin-damping tuning factor
    pc.C_span        = 0.2;   % spanwise-flow force tuning factor
    pc.C_span_torque = 0.7;   % spanwise-flow TORQUE tuning factor (independent of C_span)
    pc.k0_spanTorque = 0.2;   % reduced-frequency roll-off of the span torque
    pc.C_Tx          = 1.0;   % Tx (chordwise-axis / roll) spin-damping tuning factor
    if nargin >= 1 && isstruct(aero)
        for f = fieldnames(pc).'
            if isfield(aero, f{1}); pc.(f{1}) = aero.(f{1}); end
        end
    end
end
