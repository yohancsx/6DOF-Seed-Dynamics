function th = defaultPaperThresholds()
% DEFAULTPAPERTHRESHOLDS  Thresholds for classifyPaperMode.
%
% Separate from defaultModeThresholds on purpose: those are calibrated for OUR
% seven-mode vocabulary on OUR seed, and must not move. These are for the
% four-mode vocabulary of Hou et al. (2025) on THEIR plate.
%
% Each number is tied to something the paper states about the modes, not to a fit
% against our own output -- the point of the comparison is to find out where we
% disagree, and thresholds tuned to our trajectories would hide exactly that.
%
% Units: spins rad/s, angles deg.

    % Continuous tumbling about the LONG axis is what defines the ST family:
    % "the plate rapidly tumbling during spiral descent", ~7 rotations per spiral
    % revolution. tumbleFrac = |net rotation| / |total rotation| about body z is
    % ~1 for genuine end-over-end tumbling and ~0 for a cone that merely wobbles.
    % The discriminator is whether the plate completes a WHOLE rotation in one
    % direction. It must be measured between direction reversals, not over the
    % window: segmented ST reverses at every turning point, so its net rotation
    % over a window is near zero even though it is plainly tumbling in between.
    % A flutterer oscillates through a limited angle and never completes a turn.
    th.tumbleTurnsST = 1.0;   % complete revolutions about the long axis, one way
    th.tumbleRateST  = 5.0;   % rad/s: below this nothing is really tumbling

    % Segmented ST differs from continuous ST by REVERSALS of the tumbling
    % direction at the turning points, which show up as a much larger spread in
    % the cone angle: they measure +/-13.0 deg segmented against +/-2.3 deg
    % continuous. The split is placed between those two published values.
    th.tiltStdSegmented = 7.0;

    % Autorotation revolves about the vertical without tumbling. Their AR case
    % holds a cone of -11.4 +/- 4.2 deg, i.e. a steady attitude; the CH mode by
    % contrast has "substantial fluctuation in both velocity and orientation".
    th.vSpinAR    = 10.0;     % rad/s about the world vertical
    th.tiltStdAR  = 10.0;     % deg: cone spread above which it is not steady

    % Chaotic: "lacking consistent motion patterns", "irregular shifts between
    % phases". Detected as a run that never settles AND has no coherent rotation
    % of either kind.
    th.tiltStdChaos = 25.0;   % deg
end
