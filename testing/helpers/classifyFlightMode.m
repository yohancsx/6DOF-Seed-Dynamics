function [mode, info] = classifyFlightMode(m, th)
% CLASSIFYFLIGHTMODE  Heuristic rule-tree classifier: trajectory metrics -> mode.
%
% Vocabulary (matching the seven named working modes):
%   'gliding'      no spin, high glide ratio, low cone
%   'diving'       no spin, edge-on (high cone), fast
%   'parachuting'  no spin, broadside (low cone), slow
%   'fluttering'   flips about the span axis, no vertical circling
%   'spiral'       vertical spin + flipping, WIDE helix
%   'tightSpiral'  vertical spin + flipping, TIGHT helix
%   'autorotation' vertical spin at a STEADY cone (tiltStd ~ 0)
%   'chaotic'      did not converge
%   'undetermined' no rule matched
%
% The autorotation-vs-spiral split is the whole point: both spin about the
% vertical, but autorotation holds a steady cone (small tiltStd) while a spiral
% flips end-over-end as it circles (large tiltStd). Transparent, order-dependent
% rule tree; first match wins. Thresholds are calibrated (see defaultModeThresholds).
%
% INPUTS
%   m  : metrics struct from computeTrajectoryMetrics.
%   th : (optional) thresholds struct; defaults if omitted/empty.
%
% OUTPUTS
%   mode : char label from the vocabulary above.
%   info : struct with .reason (which rule fired), .metrics, .thresholds.

    if nargin < 2 || isempty(th); th = defaultModeThresholds(); end

    spinning = (m.verticalSpinMag > th.spinLo) || (m.spanwiseSpin > th.spinLo);
    tight    = m.helixValid && (m.helixRadius < th.helixTight);

    if ~m.converged
        mode = 'chaotic';
        reason = 'did not settle (descent/spin not converged)';

    elseif ~spinning
        % --- No rotation: glide / dive / parachute -------------------------
        if m.glideRatio > th.glideHi
            mode = 'gliding';
            reason = 'no spin, high glide ratio';
        elseif m.coneAngleDeg > th.coneEdge
            mode = 'diving';
            reason = 'no spin, edge-on (high cone)';
        elseif m.coneAngleDeg < th.coneBroad
            mode = 'parachuting';
            reason = 'no spin, broadside (low cone), slow';
        else
            mode = 'undetermined';
            reason = 'no spin, but not clearly glide/dive/parachute';
        end

    elseif (m.verticalSpinMag > th.vSpinAuto) && (m.tiltStd < th.tiltSteady)
        % --- Vertical spin at a STEADY cone -> autorotation ----------------
        mode = 'autorotation';
        reason = 'vertical-axis spin at a steady cone (low tiltStd)';

    elseif m.verticalSpinMag > th.vSpinAuto
        % --- Vertical spin + flipping cone -> spiralling -------------------
        if tight
            mode = 'tightSpiral';
            reason = 'vertical spin + flipping, tight helix';
        else
            mode = 'spiral';
            reason = 'vertical spin + flipping, wide helix';
        end

    else
        % --- Rotation without vertical circling -> flip about span --------
        mode = 'fluttering';
        reason = 'spanwise flipping, little/no vertical circling';
    end

    info.reason     = reason;
    info.metrics    = m;
    info.thresholds = th;
end
