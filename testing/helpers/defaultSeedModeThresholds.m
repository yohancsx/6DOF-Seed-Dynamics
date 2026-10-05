function th = defaultSeedModeThresholds()
% DEFAULTSEEDMODETHRESHOLDS  Every boundary between modes in classifySeedMode.
%
% These ARE the mode definitions -- testing/classifier/MODE_DEFINITIONS.md
% documents them in words and must be kept in step with this file. Each value
% says where it came from: a decision taken on the by-eye review set
% (Outputs/Classifier Review/2026-09-27_165346_3727986, decided 2026-10-03), or a
% judgement still to be checked against labels.
%
% Units: rad/s, deg, turns, spans, dimensionless ratios.

    % --- A. Rotation about the span axis (flip), from the attitude angle psi --
    th.flipNoneTurns    = 0.10;  % psi peak-to-peak below this (36 deg) -> no flip.
                                 % DECIDED 2026-10-05 (was 0.05 = 18 deg, judged):
                                 % the Hou-window revolvers rock 21-30 deg p-p
                                 % about the span while revolving tightly at a
                                 % steady tilt -- autorotation with a pitch rock,
                                 % not a flutter spiral.
    th.tumbleTurns      = 1.0;   % a one-way leg of >= 1 whole turn -> tumbling;
                                 % below it, the plate rocks -> flutter
    th.flipHysteresisTurns = 0.25; % psi must retreat this far to count a reversal
    th.periodicCV       = 0.35;  % reversal-interval CV below this -> periodic
                                 % (segmented tumbling), above -> irregular
                                 % (judgement)

    % --- B. Revolution about the vertical ------------------------------------
    th.revolveSpin      = 5.0;   % |rate of turn of the span axis's heading| above
                                 % this -> revolving (rad/s). From attitude: a
                                 % tumble or flutter about a tilted span projects
                                 % onto the body's vertical spin, but leaves the
                                 % span heading fixed. (judgement)
    th.tightRadiusSpans = 1.0;   % revolution radius < 1 span -> tight
                                 % (autorotation family); >= 1 -> wide (spiral
                                 % family). DECIDED 2026-10-03.

    % --- C. Span-axis attitude ------------------------------------------------
    th.flatSpanDeg      = 30;    % mean |theta| below this -> flat (judgement)
    th.endOnSpanDeg     = 70;    % above this -> end-on, NON-PHYSICAL.
                                 % DECIDED 2026-10-03.

    % --- D. Path, for a seed that is not revolving ----------------------------
    th.glideRatio       = 0.25;  % glide ratio above this -> glide / flutterGlide,
                                 % at or below -> parachute / flutter.
                                 % DECIDED 2026-10-03.
    th.edgeOnNormalDeg  = 45;    % plate normal more than this from vertical, no
                                 % rotation -> a dive (edge leading)
    th.edgeSlideAlong   = 0.5;   % no rotation, edge-on, span FLAT, |shat.vhat|
                                 % above this -> edge-slide, NON-PHYSICAL (draft;
                                 % none in the review set). The flat-span
                                 % condition matters: a steep span falling
                                 % straight down also lines up with v.

    % --- E. Regularity --------------------------------------------------------
    th.steadySpanStdDeg = 7.0;   % std of theta at or below this -> steady; above
                                 % -> wobbling (or chaotic, if the flip is
                                 % irregular). Same value as the paper
                                 % classifier's segmented-ST spread.
    th.chaoticSpanStdDeg = 25;   % std of theta above this -> chaotic whatever the
                                 % flip pattern (the paper classifier's CH value).
                                 % Chaotic is NON-PHYSICAL: DECIDED 2026-10-03.
end
