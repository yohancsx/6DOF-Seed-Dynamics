%% Mode ablation -- run the 7 named modes under several physics configs
% For each physics configuration, drops the seed at each mode's working nut
% position / initial condition, computes trajectory metrics, and writes a
% table to Mode_Ablation_Results.txt. Purpose: find the SIMPLEST physics config
% (fewest additions beyond the pure strip model) that still produces all modes.
%
% Read the metrics (not the classifier guess -- it's uncalibrated) to judge the
% mode: descent speed, vertical-axis spin, spanwise (tumble) spin, cone angle,
% glide ratio, helix radius.

helpersFolder = "C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\6DOF Seed Dynamics\testing\helpers";
addpath(helpersFolder);
% NOTE: physics/ and visualization/ assumed on the MATLAB path.

% --- Base config (geometry, environment, analysis) ------------------------
cfg.spanLength  = 0.050;   cfg.chordLength = 0.015;   cfg.thickness = 0.002;
cfg.bulkDensity = 65;      cfg.numStrips   = 10;      cfg.tSamples  = 0;
cfg.nutMass     = 75e-6;
cfg.rhoFluid = 1.225;   cfg.g = 9.81;
cfg.tspan = [0 10];   cfg.odeRelTol = 1e-6;   cfg.odeAbsTol = 1e-8;
cfg.metricOpts.windowStartFrac = 0.5;   cfg.metricOpts.convergeTol = 0.20;
cfg.modeThresholds = defaultModeThresholds();

xh = cfg.spanLength / 2;   yh = cfg.chordLength / 2;
baseBsp.seedShape     = polyshape([-xh, xh, xh, -xh], [-yh, -yh, yh, yh]);
baseBsp.seedDensity   = cfg.bulkDensity * cfg.thickness;
baseBsp.seedThickness = cfg.thickness;
baseBsp.numStrips     = cfg.numStrips;

c = cfg.chordLength;   S = cfg.spanLength;
qL = [1;0;0;0];   qT = axisAngleToQuat([0;0;1], pi/6);   ns = [0;0;0];

% --- The 7 modes (name, nut position, initial quaternion, initial spin) ----
modes = struct( ...
  'name', {'Spanwise-axis fluttering','Gliding','Diving','Fluttering + spiral', ...
           'Fluttering + tight spiral','Autorotation','Parachute'}, ...
  'nutPos', {[0;0;0],[0.5*c;0;0],[1.5*c;0;0],[0;0;0.01*S],[0;0;S],[c;0;1.2*S],[0;-c;0]}, ...
  'q0',     {qT, qL, qL, qT, qT, qL, qL}, ...
  'omega0', {ns, ns, ns, ns, ns, ns, ns});

% --- Physics configurations to compare ------------------------------------
% Fields: label, spanForce(SF), spanTorque(ST), geomVel(GV), copMig(CM),
%         atten(AT), Tx, C_span(Cs), C_span_torque(Ct).
% Focus: does keeping the span FORCE but dropping the span TORQUE work?
configs = struct( ...
  'label', {'FULL (working)','span torque OFF (Tx on)','span torque OFF + Tx off', ...
            'no span force (ref)','minimal strip (no span, no Tx)'}, ...
  'SF', {true, true,  true,  false, false}, ...
  'ST', {true, false, false, true,  true }, ...
  'GV', {true, true,  true,  true,  true }, ...
  'CM', {true, true,  true,  true,  true }, ...
  'AT', {false,false, false, false, false}, ...
  'Tx', {true, true,  false, true,  false}, ...
  'Cs', {0.2,  0.2,   0.2,   0.2,   0.2  }, ...
  'Ct', {0.7,  0.7,   0.7,   0.7,   0.7  });

% --- Output file ----------------------------------------------------------
outTxt = fullfile(char(fileparts(char(helpersFolder))), 'Mode_Ablation_Results.txt');
fid = fopen(outTxt, 'w');
fprintf(fid, 'Mode ablation results  %s\n', datestr(now));
fprintf(fid, ['Base seed: nut 75e-6 kg, density 65 kg/m^3, span 0.050 m, chord 0.015 m; ' ...
              'air, 10 s, from rest unless a pi/6 tilt about z is noted.\n']);
fprintf(fid, 'Judge the mode from the METRICS (classifier guess is uncalibrated).\n');
fprintf(fid, 'descent m/s | vSpin=vertical-axis spin | spanSpin=body-z (tumble) | cone deg | glide ratio | helixR m\n');

for ci = 1:numel(configs)
    cf = configs(ci);
    cfg.enableSpanForce             = cf.SF;
    cfg.enableSpanTorque            = cf.ST;
    cfg.enableSpanGeomVelocity      = cf.GV;
    cfg.enableSpanCOPMigration      = cf.CM;
    cfg.enableSpanTorqueAttenuation = cf.AT;
    cfg.enableTxDamping             = cf.Tx;
    cfg.aero = struct('C_span', cf.Cs, 'C_span_torque', cf.Ct);

    hdr = sprintf('[SF=%d ST=%d GV=%d CM=%d AT=%d Tx=%d | C_span=%.2g C_tau=%.2g]', ...
                  cf.SF, cf.ST, cf.GV, cf.CM, cf.AT, cf.Tx, cf.Cs, cf.Ct);
    fprintf(fid, '\n===========================================================================================\n');
    fprintf(fid, 'CONFIG: %s   %s\n', cf.label, hdr);
    fprintf(fid, '-------------------------------------------------------------------------------------------\n');
    fprintf(fid, '%-26s %-13s %7s %7s %8s %5s %6s %8s %4s\n', ...
            'mode', 'guess', 'descent', 'vSpin', 'spanSpin', 'cone', 'glide', 'helixR', 'conv');
    fprintf('\n=== %s %s ===\n', cf.label, hdr);

    for mi = 1:numel(modes)
        md = modes(mi);
        try
            r = runSingleMode(md.name, md.nutPos, md.q0, md.omega0, cfg, baseBsp);
            m = r.metrics;
            hR = m.helixRadius; if ~m.helixValid; hR = NaN; end
            fprintf(fid, '%-26s %-13s %7.2f %7.2f %8.2f %5.0f %6.2f %8.3f %4d\n', ...
                    md.name, r.mode, m.descentSpeed, m.verticalSpinMag, m.spanwiseSpin, ...
                    m.coneAngleDeg, m.glideRatio, hR, m.converged);
        catch ME
            fprintf(fid, '%-26s FAILED: %s\n', md.name, ME.message);
            fprintf('  %s FAILED: %s\n', md.name, ME.message);
        end
    end
end

fclose(fid);
fprintf('\nWrote results to: %s\n', outTxt);
