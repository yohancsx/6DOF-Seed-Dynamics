function th = defaultModeThresholds()
% DEFAULTMODETHRESHOLDS  Tunable thresholds for classifyFlightMode.
%
% These are STARTING GUESSES, not calibrated values. The intended workflow:
% run runSeedModeSuite on cases whose mode you already know, compare the
% classifier output to reality, and adjust these until the labels match --
% THEN freeze them. Get a copy to edit with:
%     th = defaultModeThresholds();  th.glideHi = 0.7;  % etc.
% and pass it to classifyFlightMode / runModeSweep.
%
% Units: spins in rad/s, cone/tilt in degrees, speeds in m/s, radii in m,
% glide ratio and straightness dimensionless.

    th.spinVhi    = 5.0;   % "fast" vertical-axis spin -> autorotation
    th.spinVlo    = 1.0;   % "some" vertical-axis spin -> spiral
    th.tumbleHi   = 5.0;   % sustained end-over-end (spanwise) spin -> tumbling
    th.tumbleLo   = 1.0;   % below this, end-over-end is negligible
    th.coneLo     = 20;    % deg: plate near-flat / broadside
    th.coneHi     = 60;    % deg: plate near edge-on
    th.glideHi    = 0.5;   % glide ratio above which motion counts as gliding
    th.glideLo    = 0.2;   % glide ratio below which it is not gliding
    th.straightHi = 0.8;   % straightness (0..1) required for gliding
    th.descentLo  = 1.0;   % "slow" descent (autorotation / parachute)
    th.descentHi  = 2.0;   % "fast" descent (dive)
    th.helixTight = 0.05;  % m: helix radius below which the spiral is "tight"
    th.flutterTilt= 10;    % deg: tilt std above which the plate is oscillating
end
