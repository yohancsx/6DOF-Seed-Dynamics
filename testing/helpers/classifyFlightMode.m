function [mode, info] = classifyFlightMode(m, th)
% CLASSIFYFLIGHTMODE  Heuristic rule-tree classifier: trajectory metrics -> mode.
%
% Maps the descriptors from computeTrajectoryMetrics to a named flight mode.
% Vocabulary (informed by the samara / falling-plate literature):
%   'gliding'        sloped, constant-speed, low spin, straight
%   'diving'         steep, edge-on, fast
%   'parachuting'    broadside, slow, little horizontal travel (CoM-below)
%   'tumbling'       sustained end-over-end (spanwise-axis) rotation
%   'spiralTumbling' vertical-axis spin along a WIDE helix (large horizontal spread)
%   'tightSpiral'    vertical-axis spin along a TIGHT helix
%   'autorotation'   steady vertical-axis spin, small cone, slow, near-centred
%   'fluttering'     oscillating tilt, no net spin (side-to-side falling card)
%   'chaotic'        did not settle / irregular
%   'undetermined'   no rule matched (inspect the metrics)
%
% This is a TRANSPARENT, order-dependent rule tree -- the first matching branch
% wins. Thresholds are tunable (see defaultModeThresholds); calibrate them
% against known cases before trusting the output.
%
% INPUTS
%   m  : metrics struct from computeTrajectoryMetrics.
%   th : (optional) thresholds struct; defaults if omitted/empty.
%
% OUTPUTS
%   mode : char label from the vocabulary above.
%   info : struct with .reason (which rule fired), .metrics, .thresholds.

    if nargin < 2 || isempty(th); th = defaultModeThresholds(); end

    spinV = m.verticalSpinMag;   % vertical-axis spin  (spiral / autorotation)
    spinZ = m.spanwiseSpin;      % end-over-end spin   (tumbling)
    cone  = m.coneAngleDeg;      % 0 = flat/broadside, 90 = edge-on
    gr    = m.glideRatio;
    tight = m.helixValid && m.helixRadius < th.helixTight;

    if ~m.converged
        mode = 'chaotic';
        reason = 'did not settle (descent/spin not converged)';

    elseif spinV > th.spinVhi && cone < th.coneLo ...
            && m.descentSpeed < th.descentLo && tight
        mode = 'autorotation';
        reason = 'steady fast vertical spin, small cone, slow, near-centred';

    elseif spinZ > th.tumbleHi
        mode = 'tumbling';
        reason = 'sustained end-over-end (spanwise) rotation';

    elseif spinV > th.spinVlo && m.helixValid
        if tight
            mode = 'tightSpiral';
            reason = 'vertical-axis spin along a tight helix';
        else
            mode = 'spiralTumbling';
            reason = 'vertical-axis spin along a wide helix';
        end

    elseif gr > th.glideHi && spinV < th.spinVlo ...
            && spinZ < th.tumbleLo && m.straightness > th.straightHi
        mode = 'gliding';
        reason = 'high glide ratio, straight track, low spin';

    elseif gr < th.glideLo && cone > th.coneHi && m.descentSpeed > th.descentHi
        mode = 'diving';
        reason = 'steep, edge-on, fast descent';

    elseif gr < th.glideLo && cone < th.coneLo && m.descentSpeed < th.descentLo
        mode = 'parachuting';
        reason = 'broadside, slow, little horizontal travel';

    elseif m.tiltStd > th.flutterTilt && spinV < th.spinVlo && spinZ < th.tumbleLo
        mode = 'fluttering';
        reason = 'oscillating tilt with no net spin';

    else
        mode = 'undetermined';
        reason = 'no rule matched -- inspect metrics';
    end

    info.reason     = reason;
    info.metrics    = m;
    info.thresholds = th;
end
