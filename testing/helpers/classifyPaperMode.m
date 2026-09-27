function [mode, info] = classifyPaperMode(pm, th)
% CLASSIFYPAPERMODE  Trajectory metrics -> the four modes of Hou et al. (2025).
%
% Vocabulary (their Fig. 2a):
%   'AR'   Autorotation    -- stable descent revolving about a vertical axis, with
%                             NO tumbling and a steady, shallow cone.
%   'CST'  Continuous ST   -- continuous rapid tumbling about the plate's LONGER
%                             axis while spiralling down a vertical axis, at a
%                             nearly constant cone angle.
%   'SST'  Segmented ST    -- the same, but with distinct turning points where the
%                             tumbling direction REVERSES, and a much larger cone
%                             spread (+/-13.0 deg vs +/-2.3 deg).
%   'CH'   Chaotic         -- irregular tumbling and oscillation, no sustained
%                             periodic pattern, large fluctuation in orientation.
%   'FA'   Falling         -- rapid, predominantly vertical descent at a large
%                             tilt, heavier side down, without sustained rotation.
%   'FL'   Fluttering      -- NOT one of their modes. Rocking fast about the long
%                             axis without ever completing a turn, at a flat
%                             attitude. Their plate never does this; our model
%                             does, over much of the window, so it gets its own
%                             label rather than being folded into FA.
%
% This is DELIBERATELY a separate classifier from classifyFlightMode: that one owns
% our own seven-mode vocabulary and its thresholds are calibrated for our seed, and
% neither should move to accommodate the other. Report both side by side.
%
% THE DISCRIMINATORS, and why these and not others:
%   selfTumbleFrac  separates genuine end-over-end tumbling about the long axis
%                 (~1) from a cone that merely wobbles (~0). This is the ST/AR
%                 split, and it is a ratio, so it needs no velocity scale. It is
%                 built on psi' -- omega_z with the REVOLUTION's projection on the
%                 tilted span axis removed. Using omega_z raw would label a pure
%                 revolution at a tilt as tumbling (testPaperMetrics case A).
%   tumbleReversals + spanAxisTiltStd   separate segmented from continuous ST --
%                 their turning points are literally reversals of the tumbling
%                 direction, and the paper's own cone spreads differ 5-fold.
%   verticalSpinMag   detects revolution about the vertical (AR).
%   spanAxisTiltStd   detects the "substantial fluctuation in orientation" that
%                 defines CH, and the steady attitude that defines AR.
% Descent SPEED is deliberately NOT used: the paper reports FA as fastest and AR
% as slowest, so using speed to define them would make that a tautology instead of
% a prediction we can check.
%
% INPUTS
%   pm : metrics struct from computePaperMetrics (needs .selfTumbleFrac,
%        .selfRotationRate, .tumbleReversals, .spanAxisTiltStd, .verticalSpinMag,
%        .converged).
%   th : (optional) thresholds; defaultPaperThresholds() if omitted.
%
% OUTPUTS
%   mode : 'AR' | 'CST' | 'SST' | 'CH' | 'FA' | 'FL'
%   info : struct with .reason (which rule fired), .metrics, .thresholds.

    if nargin < 2 || isempty(th); th = defaultPaperThresholds(); end

    tumbling  = (pm.maxTurnsOneWay > th.tumbleTurnsST) && (pm.selfRotationRate > th.tumbleRateST);
    revolving =  pm.verticalSpinMag > th.vSpinAR;
    steady    =  pm.spanAxisTiltStd < th.tiltStdAR;

    if tumbling
        % --- The ST family: whole rotations about the long axis
        if pm.tumbleReversals > 0 || pm.spanAxisTiltStd > th.tiltStdSegmented
            mode   = 'SST';
            reason = 'tumbling about the long axis, with direction reversals / large cone spread';
        else
            mode   = 'CST';
            reason = 'continuous tumbling about the long axis at a steady cone';
        end

    elseif pm.spanAxisTiltStd > th.tiltStdChaos
        % --- The attitude never settles. Their CH has "substantial fluctuation in
        % both velocity and orientation"; FA by contrast holds its large tilt
        % CONSISTENTLY. So a wildly swinging attitude is CH whether or not the
        % descent speed happens to have plateaued.
        mode   = 'CH';
        reason = 'orientation never settles (large cone swing), no sustained pattern';

    elseif ~pm.converged && pm.spanAxisTiltStd > th.tiltStdSegmented
        mode   = 'CH';
        reason = 'unsettled, with an unsteady cone and no sustained pattern';

    elseif revolving && steady
        % --- Revolving about the vertical at a steady attitude, without tumbling
        mode   = 'AR';
        reason = 'vertical-axis revolution at a steady cone, no tumbling';

    elseif revolving
        % --- Rotating, but the attitude never settles
        mode   = 'CH';
        reason = 'rotating but the cone never settles';

    elseif pm.selfRotationRate > th.tumbleRateST && pm.maxTurnsOneWay < th.tumbleTurnsST
        % --- FL: rocking hard about the long axis without ever completing a turn.
        % DELIBERATELY OUTSIDE their four-mode vocabulary. Their plate never does
        % this, so there is no published label for it -- and forcing it into FA
        % would be wrong twice over: FA is defined by a LARGE, CONSISTENT tilt and
        % no sustained rotation, whereas this is a flat plate rotating fast enough
        % to reverse direction hundreds of times. Reporting it as its own label
        % keeps the disagreement with the paper visible instead of hiding it in a
        % catch-all.
        mode   = 'FL';
        reason = 'fluttering: fast rocking about the long axis, never a full turn';

    else
        % --- Descending without sustained rotation of either kind
        mode   = 'FA';
        reason = 'no tumbling and no coherent revolution: a tilted fall';
    end

    info.reason     = reason;
    info.metrics    = pm;
    info.thresholds = th;
end
