function pal = seedModePalette()
% SEEDMODEPALETTE  Names, one-letter codes, colours and physicality of the
%   six-part classifier's modes (classifySeedMode), for maps and ASCII grids.
%
% Order = the precedence order of testing/classifier/MODE_DEFINITIONS.md section
% 3, plus 'failed' (integration error or non-finite state) last. Non-physical
% modes are drawn in reds/greys so they stand out on a map; physical families
% share a hue (dives brown, glides green, revolvers blue, flutters orange,
% tumbles purple).
%
% OUTPUT pal (struct):
%   .names    1x18 cellstr
%   .codes    1x18 char, one letter per mode (ASCII maps)
%   .colors   18x3 RGB in [0,1]
%   .physical 1x18 logical ('failed' is false)
%   .paper    1x18 cellstr, the Hou et al. (2025) family each mode corresponds
%             to: 'AR' | 'ST' | 'CH' | 'FA' | '' (no equivalent), as in the
%             "paper" column of MODE_DEFINITIONS section 3. NOTE chaotic -> CH is
%             recorded for information, but summarizeModeGrid still scores every
%             non-physical label as a DISAGREEMENT (the documented rule) and only
%             reports the chaotic-in-CH count separately.

    T = { ...
    % name                 code  RGB                      physical  paper
      'endOnFall',          'e', [0.55 0.00 0.00],        false,    ''
      'endOnSpin',          'E', [0.85 0.10 0.10],        false,    ''
      'chaotic',            'X', [0.25 0.25 0.25],        false,    'CH'
      'wobblingSpin',       'W', [0.95 0.45 0.45],        false,    ''
      'edgeSlide',          'Z', [0.60 0.60 0.60],        false,    ''
      'chordDive',          'd', [0.55 0.35 0.15],        true,     ''
      'obliqueDive',        'D', [0.80 0.55 0.25],        true,     'FA'
      'glide',              'g', [0.40 0.75 0.35],        true,     ''
      'parachute',          'p', [0.70 0.90 0.55],        true,     ''
      'autorotation',       'A', [0.10 0.30 0.85],        true,     'AR'
      'spiralGlide',        'G', [0.45 0.70 0.95],        true,     ''
      'flutter',            'f', [1.00 0.70 0.20],        true,     ''
      'flutterGlide',       'F', [1.00 0.85 0.45],        true,     ''
      'flutterSpiral',      'L', [0.95 0.55 0.05],        true,     ''
      'tumble',             't', [0.75 0.60 0.90],        true,     ''
      'tightSpiralTumble',  'T', [0.45 0.20 0.70],        true,     'ST'
      'spiralTumble',       'S', [0.65 0.40 0.85],        true,     'ST'
      'failed',             '?', [1.00 1.00 1.00],        false,    ''
    };
    pal.names    = T(:,1).';
    pal.codes    = [T{:,2}];
    pal.colors   = vertcat(T{:,3});
    pal.physical = [T{:,4}];
    pal.paper    = T(:,5).';
end
