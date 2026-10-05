function [mode, parts, info] = classifySeedMode(sm, th)
% CLASSIFYSEEDMODE  Six-part flight-mode classifier (the "new" classifier).
%
% Judges six things INDEPENDENTLY, then builds the mode name from them, so two
% runs that move the same way always get the same name. The parts and every
% boundary are defined in testing/classifier/MODE_DEFINITIONS.md, and the
% numbers live in defaultSeedModeThresholds -- keep all three in step.
%
%   A  flip        rotation about the span axis:
%                  'none' | 'flutter' | 'tumble' | 'segmentedTumble' | 'irregular'
%   B  revolution  about the world vertical:  'none' | 'tight' | 'wide'
%   C  attitude    span axis from horizontal: 'flat' | 'steep' | 'endOn'
%   D  path        non-revolving seeds only:  'steep' | 'glide' | '' (n/a)
%   E  regularity  'steady' | 'wobbling' | 'chaotic'
%   F  physical    true | false
%
% MODE NAMES, in order of precedence (first match wins):
%   endOnFall    C endOn, no rotation                         NON-PHYSICAL
%   endOnSpin    C endOn, flipping or revolving               NON-PHYSICAL
%   chaotic      E chaotic (irregular whole-turn reversals, or
%                theta std > th.chaoticSpanStdDeg)            NON-PHYSICAL
%   wobblingSpin E wobbling, B revolving                      NON-PHYSICAL
%   edgeSlide    no rotation, edge-on, C flat, moving along
%                the span                                     NON-PHYSICAL
%   chordDive    no rotation, edge-on, C flat
%   obliqueDive  no rotation, edge-on, C steep
%   glide        no rotation, not edge-on, glide ratio > th.glideRatio
%   parachute    no rotation, not edge-on, glide ratio <= th.glideRatio
%   autorotation A none,    B tight
%   spiralGlide  A none,    B wide
%   flutter      A flutter, B none, glide ratio <= th.glideRatio
%   flutterGlide A flutter, B none, glide ratio > th.glideRatio
%   flutterSpiral A flutter, B tight or wide
%   tumble       A tumble/segmentedTumble, B none
%   tightSpiralTumble  A tumble/segmentedTumble, B tight
%   spiralTumble       A tumble/segmentedTumble, B wide
%
% INPUTS
%   sm : struct from computeSeedModeMetrics.
%   th : (optional) thresholds; defaultSeedModeThresholds() if omitted.
%
% OUTPUTS
%   mode  : char, one of the names above.
%   parts : struct .flip .revolution .attitude .path .regularity .physical
%   info  : struct .reason, .metrics, .thresholds

    if nargin < 2 || isempty(th); th = defaultSeedModeThresholds(); end

    % --- A. flip -----------------------------------------------------------
    if sm.flipPPTurns < th.flipNoneTurns
        A = 'none';
    elseif sm.flipMaxTurnsOneWay < th.tumbleTurns
        A = 'flutter';
    elseif sm.flipReversals == 0
        A = 'tumble';
    elseif sm.flipReversals < 3 || sm.flipReversalCV < th.periodicCV
        A = 'segmentedTumble';
    else
        A = 'irregular';
    end

    % --- B. revolution -----------------------------------------------------
    if abs(sm.spanHeadingRate) <= th.revolveSpin
        B = 'none';
    elseif ~isfinite(sm.radiusPerSpan) || sm.radiusPerSpan < th.tightRadiusSpans
        B = 'tight';          % no valid helix while revolving = turning in place
    else
        B = 'wide';
    end

    % --- C. attitude -------------------------------------------------------
    if sm.absSpanTiltDeg > th.endOnSpanDeg
        C = 'endOn';
    elseif sm.absSpanTiltDeg < th.flatSpanDeg
        C = 'flat';
    else
        C = 'steep';
    end

    % --- D. path (non-revolving only) --------------------------------------
    if strcmp(B, 'none')
        if sm.glideRatio > th.glideRatio; D = 'glide'; else; D = 'steep'; end
    else
        D = '';
    end

    % --- E. regularity -----------------------------------------------------
    if strcmp(A, 'irregular') || sm.spanAxisTiltStd > th.chaoticSpanStdDeg
        E = 'chaotic';
    elseif sm.spanAxisTiltStd > th.steadySpanStdDeg
        E = 'wobbling';
    else
        E = 'steady';
    end

    rotating = ~strcmp(A, 'none') || ~strcmp(B, 'none');
    edgeOn   = sm.coneAngleDeg > th.edgeOnNormalDeg;
    tumbling = any(strcmp(A, {'tumble', 'segmentedTumble'}));

    % --- name, in order of precedence --------------------------------------
    if strcmp(C, 'endOn') && ~rotating
        mode = 'endOnFall';     reason = 'span axis near vertical, no rotation';
    elseif strcmp(C, 'endOn')
        mode = 'endOnSpin';     reason = 'span axis near vertical, rotating';
    elseif strcmp(E, 'chaotic')
        mode = 'chaotic';       reason = 'whole turns about the span with irregular reversals';
    elseif strcmp(E, 'wobbling') && ~strcmp(B, 'none')
        mode = 'wobblingSpin';  reason = 'revolving, but the span axis wobbles';
    elseif ~rotating && edgeOn && strcmp(C, 'flat') && sm.spanAlongVel > th.edgeSlideAlong
        mode = 'edgeSlide';     reason = 'no rotation, edge-on, sliding along its level span';
    elseif ~rotating && edgeOn && strcmp(C, 'flat')
        mode = 'chordDive';     reason = 'no rotation, edge-on, span level';
    elseif ~rotating && edgeOn
        mode = 'obliqueDive';   reason = 'no rotation, edge-on, span tilted';
    elseif ~rotating && strcmp(D, 'glide')
        mode = 'glide';         reason = 'no rotation, travels sideways';
    elseif ~rotating
        mode = 'parachute';     reason = 'no rotation, falls nearly straight down';
    elseif strcmp(A, 'none') && strcmp(B, 'tight')
        mode = 'autorotation';  reason = 'revolves within a span of the seed, no flip';
    elseif strcmp(A, 'none')
        mode = 'spiralGlide';   reason = 'revolves on a radius of a span or more, no flip';
    elseif strcmp(A, 'flutter') && strcmp(B, 'none') && strcmp(D, 'glide')
        mode = 'flutterGlide';  reason = 'rocks about the span, travels sideways';
    elseif strcmp(A, 'flutter') && strcmp(B, 'none')
        mode = 'flutter';       reason = 'rocks about the span, falls nearly straight down';
    elseif strcmp(A, 'flutter')
        mode = 'flutterSpiral'; reason = 'rocks about the span while revolving';
    elseif tumbling && strcmp(B, 'none')
        mode = 'tumble';        reason = 'whole turns about the span, no revolution';
    elseif tumbling && strcmp(B, 'tight')
        mode = 'tightSpiralTumble'; reason = 'tumbles while revolving within a span';
    else
        mode = 'spiralTumble';  reason = 'tumbles while revolving on a radius of a span or more';
    end

    physical = ~any(strcmp(mode, {'endOnFall', 'endOnSpin', 'chaotic', ...
                                  'wobblingSpin', 'edgeSlide'}));
    parts = struct('flip', A, 'revolution', B, 'attitude', C, 'path', D, ...
                   'regularity', E, 'physical', physical);
    info  = struct('reason', reason, 'metrics', sm, 'thresholds', th);
end
