function th = defaultModeThresholds()
% DEFAULTMODETHRESHOLDS  Thresholds for classifyFlightMode.
%
% Calibrated against the seven known working modes (see the mode-ablation /
% metric-dump runs) with the default "FULL-minus-geomVelocity" physics. The key
% discriminator is tiltStd (cone-angle variation): a spinning seed with a STEADY
% cone (tiltStd ~ 0) is autorotating; one whose cone swings (tiltStd large) is
% tumbling/spiralling. Re-calibrate if the physics or seed changes.
%
% Units: spins rad/s, angles/cone deg, radii m, glide ratio dimensionless.

    th.spinLo     = 5.0;   % below this on BOTH spin axes -> not rotating
    th.vSpinAuto  = 10;    % vertical-axis spin needed for autorotation / spiral
    th.tiltSteady = 5.0;   % deg: cone-angle std below which the cone is "steady"
                           % (autorotation) vs. flipping (tumbling/spiral)
    th.tiltChaos  = 30;    % deg: cone-angle std ABOVE which an unsettled, non-
                           % spinning run counts as chaotic (genuine incoherent
                           % tumbling) rather than a slow-settling named mode.
                           % Gates 'chaotic' so non-convergence alone no longer
                           % dumps coherent spirals/dives into 'chaotic'.
    th.helixTight = 0.05;  % m: horizontal helix radius below which a spiral is "tight"
    th.glideHi    = 1.0;   % glide ratio above which motion counts as gliding
    th.coneEdge   = 45;    % deg: cone above this (no spin) -> diving (edge-on)
    th.coneBroad  = 20;    % deg: cone below this (no spin) -> parachuting (broadside)
end
