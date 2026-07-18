function [cases, baselineBsp] = seedTestCases(cfg)
% SEEDTESTCASES  Build the full list of test-suite cases (as data) plus the
%   baseline base-seed-params they are derived from.
%
% Each case is a struct describing ONE simulation:
%   .group        - group name (subfolder + overlay grouping), e.g. 'nutMass'
%   .label        - filesystem-safe case name (figure title + filename)
%   .bspOverrides - fields to override in the baseline base-seed-params before
%                   rebuilding via setupSeedShapeAndMass, or [] to reuse the
%                   cached baseline seed (IC-only cases)
%   .postBuildFn  - @(sp)->sp mutation applied AFTER the seed is built (used for
%                   per-strip lift/drag multipliers), or []
%   .q0,.v0,.omega0 - initial orientation (quat), inertial velocity, body
%                   angular velocity. Default: dropped from rest, level.
%   .sweepVal     - the swept scalar (for reference/ordering)
%
% Body-axis convention (x=chord, y=normal/vertical, z=span), matching the model:
%   pitch = rotation about body z (spanwise); roll = about body x (chordwise);
%   yaw spin = angular velocity about body y (normal) -- the autorotation axis.
%
% INPUT
%   cfg : suite config (geometry, sweep magnitudes, increments, group toggles).
%
% OUTPUTS
%   cases       : 1xN struct array of case structs (only enabled groups).
%   baselineBsp : the baseline base-seed-params struct.

% -------------------------------------------------------------------------
% Baseline base-seed-params (the symmetric rectangular reference seed).
% -------------------------------------------------------------------------
baselineBsp.seedShape     = rectanglePolyshape(cfg.spanLength, cfg.chordLength);
baselineBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;   % areal density (kg/m^2)
baselineBsp.seedThickness = cfg.thickness;
baselineBsp.numStrips     = cfg.numStrips;
baselineBsp.tSamples      = cfg.tSamples;
baselineBsp.nutMass_t     = cfg.nutMass * ones(size(cfg.tSamples));
baselineBsp.nutPos_t      = repmat(cfg.nutPos, 1, numel(cfg.tSamples));

% -------------------------------------------------------------------------
% Assemble enabled groups. Each generator returns a cell array of case structs;
% concatenating cells then [C{:}] yields one struct array (all cases share the
% same field set from baseCase).
% -------------------------------------------------------------------------
g = cfg.groups;
C = {};
if g.nutMass;   C = [C, genNutMass(cfg, baselineBsp)];        end
if g.nutChord;  C = [C, genNutPos(cfg, baselineBsp, 'chord')];end
if g.nutSpan;   C = [C, genNutPos(cfg, baselineBsp, 'span')]; end
if g.nutDiag;   C = [C, genNutPos(cfg, baselineBsp, 'diag')]; end
if g.pitch;     C = [C, genTilt(cfg, 'pitch')];              end
if g.roll;      C = [C, genTilt(cfg, 'roll')];               end
if g.yawSpin;   C = [C, genYawSpin(cfg)];                    end
if g.asymmetry; C = [C, genAsymmetry(cfg)];                  end
if g.stripConv; C = [C, genStripConv(cfg)];                 end

cases = [C{:}];
end


% =========================================================================
% GROUP GENERATORS (each returns a cell array of case structs)
% =========================================================================

function C = genNutMass(cfg, baselineBsp)
% Vary the nut mass +/- cfg.nutMassFrac about baseline, in cfg.nIncr steps.
    C = {};
    fracs        = linspace(-cfg.nutMassFrac, cfg.nutMassFrac, cfg.nIncr);
    baseNutMass  = baselineBsp.nutMass_t(1);
    for f = fracs
        c = baseCase();
        c.group        = 'nutMass';
        c.label        = sprintf('nutMass_%+dpct', round(100*f));
        newMass        = baseNutMass * (1 + f);
        c.bspOverrides = struct('nutMass_t', newMass * ones(size(baselineBsp.tSamples)));
        c.sweepVal     = newMass;
        C{end+1} = c; %#ok<AGROW>
    end
end

function C = genNutPos(cfg, baselineBsp, axisName)
% Move the nut from centre outward along chord / span / diagonal, from 0 to
% cfg.nutPosMaxFrac of the half-dimension (values >1 place the nut OFF the
% planform -- allowed on purpose, mirrors some edge-case samaras).
    C  = {};
    fracs = linspace(0, cfg.nutPosMaxFrac, cfg.nIncr);
    hc = cfg.chordLength / 2;   % half-chord (body x)
    hs = cfg.spanLength  / 2;   % half-span  (body z)
    for f = fracs
        c = baseCase();
        switch axisName
            case 'chord'
                pos = [f*hc; 0; 0];        c.group = 'nutChord';
                c.label = sprintf('nutChord_%03.0fpct', round(100*f));
            case 'span'
                pos = [0; 0; f*hs];        c.group = 'nutSpan';
                c.label = sprintf('nutSpan_%03.0fpct', round(100*f));
            case 'diag'
                pos = [f*hc; 0; f*hs];     c.group = 'nutDiag';
                c.label = sprintf('nutDiag_%03.0fpct', round(100*f));
        end
        c.bspOverrides = struct('nutPos_t', repmat(pos, 1, numel(baselineBsp.tSamples)));
        c.sweepVal     = f;
        C{end+1} = c; %#ok<AGROW>
    end
end

function C = genTilt(cfg, axisName)
% Drop from an initial pitch or roll tilt, -cfg.tiltMaxDeg .. +cfg.tiltMaxDeg.
% IC-only (reuse the baseline seed). +theta vs -theta should mirror for the
% symmetric seed -- the built-in symmetry check.
    C = {};
    degs = linspace(-cfg.tiltMaxDeg, cfg.tiltMaxDeg, cfg.nIncr);
    switch axisName
        case 'pitch'; ax = [0; 0; 1]; grp = 'pitch';   % about body z (spanwise)
        case 'roll';  ax = [1; 0; 0]; grp = 'roll';    % about body x (chordwise)
    end
    for d = degs
        c = baseCase();
        c.group    = grp;
        c.label    = sprintf('%s_%+03.0fdeg', grp, d);
        c.q0       = axisAngleToQuat(ax, deg2rad(d));
        c.sweepVal = d;
        C{end+1} = c; %#ok<AGROW>
    end
end

function C = genYawSpin(cfg)
% Drop with an initial spin about the plate normal (body y) -- the autorotation
% axis -- from cfg.yawSpinMin to cfg.yawSpinMax rad/s. IC-only.
    C = {};
    rates = linspace(cfg.yawSpinMin, cfg.yawSpinMax, cfg.nIncr);
    for w = rates
        c = baseCase();
        c.group    = 'yawSpin';
        c.label    = sprintf('yawSpin_%.1frads', w);
        c.omega0   = [0; w; 0];
        c.sweepVal = w;
        C{end+1} = c; %#ok<AGROW>
    end
end

function C = genAsymmetry(cfg)
% Halve lift+drag / lift-only / drag-only on one spanwise side (z>0) using the
% per-strip multipliers. Applied post-build (no rebuild) via postBuildFn.
    C = {};
    modes = {'both', 'lift', 'drag'};
    for i = 1:numel(modes)
        m = modes{i};
        c = baseCase();
        c.group       = 'asymmetry';
        c.label       = sprintf('asym_half%s', [upper(m(1)), m(2:end)]);
        c.postBuildFn = @(sp) applyAsymmetry(sp, cfg.asymFactor, m);
        c.sweepVal    = i;
        C{end+1} = c; %#ok<AGROW>
    end
end

function C = genStripConv(cfg)
% Strip-count convergence check: same drop-from-rest case at each numStrips in
% cfg.stripCounts. Trajectories should converge as the count increases.
    C = {};
    for n = cfg.stripCounts
        c = baseCase();
        c.group        = 'stripConv';
        c.label        = sprintf('strips_%02d', n);
        c.bspOverrides = struct('numStrips', n);
        c.sweepVal     = n;
        C{end+1} = c; %#ok<AGROW>
    end
end


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function c = baseCase()
% Default case: dropped from rest, level, no spin, reusing the baseline seed.
    c.group        = '';
    c.label        = '';
    c.bspOverrides = [];            % [] -> reuse cached baseline seedParams
    c.postBuildFn  = [];            % [] -> no post-build mutation
    c.q0           = [1; 0; 0; 0];  % identity orientation (flat, broadside drop)
    c.v0           = [0; 0; 0];     % from rest (inertial)
    c.omega0       = [0; 0; 0];     % no angular velocity (body)
    c.sweepVal     = 0;
end

function ps = rectanglePolyshape(spanLen, chordLen)
% Rectangular wing in the drawing frame (x = span -> body z, y = chord -> body x),
% centred at the origin. Matches the setup convention in setupSeedShapeAndMass.
    xh = spanLen  / 2;
    yh = chordLen / 2;
    ps = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
end

function sp = applyAsymmetry(sp, factor, mode)
% Scale lift and/or drag by `factor` on the z>0 spanwise half (midline split).
    sideMask = sp.strips.zgc_body > 0;   % "one side" = spanwise midline split
    if any(strcmp(mode, {'both', 'lift'}))
        sp.strips.liftMult(sideMask) = factor;
    end
    if any(strcmp(mode, {'both', 'drag'}))
        sp.strips.dragMult(sideMask) = factor;
    end
end
