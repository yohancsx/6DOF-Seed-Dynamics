function [modeList, modeCol, modeCode] = seedModeColors()
% SEEDMODECOLORS  Canonical flight-mode vocabulary + display colours.
%
% One source of truth for the mode names, their single-letter codes, and the
% RGB colour each mode is drawn in, so every mode plot/animation (the phase map,
% the mode-coloured trajectory animation, ...) renders a given mode identically.
% The row order defines the integer index used by classifyModeTimeline and the
% phase-map grid.
%
% OUTPUTS
%   modeList : 1x10 cell of char labels (matches classifyFlightMode's vocabulary
%              plus 'failed' for integrator blow-ups).
%   modeCol  : 10x3 RGB, one row per mode (same order as modeList).
%   modeCode : 1x10 char, one single-letter code per mode (ASCII summaries).

    modeList = {'gliding','diving','parachuting','fluttering','spiral', ...
                'tightSpiral','autorotation','chaotic','undetermined','failed'};
    modeCode = 'gdpfstacux';
    modeCol  = [0.20 0.70 0.20;   % gliding      green
                0.70 0.10 0.10;   % diving       dark red
                0.20 0.80 0.80;   % parachuting  cyan
                0.95 0.60 0.10;   % fluttering   orange
                0.20 0.40 0.90;   % spiral       blue
                0.55 0.20 0.70;   % tightSpiral  purple
                0.92 0.20 0.60;   % autorotation magenta
                0.30 0.30 0.30;   % chaotic      dark grey
                0.75 0.75 0.75;   % undetermined light grey
                0.00 0.00 0.00];  % failed       black
end
