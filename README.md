# 6-DOF Seed Dynamics
Created by Yohan Sequeira, with great help from Olivia Pomerenk

A MATLAB simulation of the free flight of an autorotating winged seed (a **samara**),
modelled as a thin plate carrying a discrete "nut" mass, using **quasi-steady strip-theory
aerodynamics**. It is a full 6-degree-of-freedom (13-state) extension of the 2D
Andersen–Pesavento–Wang / Pomerenk–Ristroph falling-plate model into three dimensions.

The goal is to reproduce the diverse descent behaviours real seeds show — gliding, diving,
fluttering, tumbling/spiralling, autorotation — from **where the seed's mass sits** and
**what shape its wing has**, while adding as little physics that is not independently
justified as possible.

> **How we got here:** [`RESEARCH_LOG.md`](RESEARCH_LOG.md) records what was tried, in what
> order, what it showed, and why the direction changed. This README describes the code as it
> is now.

---

## Current state (September 2026)

There are **two models** sharing one rigid-body core:

| model | builder / RHS | status |
|---|---|---|
| `planar` | `setupSeedShapeAndMass` / `seed6DOFODE` | **Frozen reference.** Flat seed. Reproduces all seven reference modes, but only with four empirical whole-seed terms (span force, span torque with migrating CoP, `Tx`, `Ty`) carrying tuned constants. Kept byte-identical so every earlier run stays comparable. |
| `shape3d` | `setupSeedShape3D` / `seed6DOFODE3D` | **Active model.** Spanwise twist and curvature; the empirical whole-seed terms removed; first-principles added mass (Kirchhoff's pair) and tip edge drag on by default; LEV lift available but off. |

What the `shape3d` model gets right and wrong at present (details and numbers in the log,
phases 5–6):

- **Right:** descent speed on the Hou et al. (2025) plate, 1.04–1.36 × `sqrt(σg/ρ)` against
  a published O(1) prefactor.
- **Wrong:** attitude. The model cannot hold a fixed attitude about its long axis, and 0/104
  cells agree with the published phase map. *Both measured before a sign error in the paper
  metrics was fixed (2026-10-03). That error made every steady revolver read as spiral
  tumbling, so the comparison has to be re-run (log, phase 7).*
- **Wrong:** with spanwise flow. Without the planar patches it degenerates by rolling and
  sliding spanwise, and can fall span-down.

The diagnosis is missing low-aspect-ratio (finite-wing) aerodynamics. The plan to supply it
is in the [roadmap](#roadmap).

---

## Conventions (worth knowing before reading the code)

- **Body axes:** `x` = chordwise, `y` = plate-normal, `z` = spanwise. A flat plate lies in the
  body `x`–`z` plane.
- **World frame:** Y-up. Gravity is `[0; -g; 0]`.
- **State vector (13×1):** `x = [ r(3) ; q(4) ; v(3) ; omega(3) ]`
  - `r`: CoM position, **inertial** frame
  - `q`: orientation quaternion `[q0;q1;q2;q3]`, scalar-first, **body → world**
  - `v`: CoM velocity, **inertial** frame
  - `omega`: angular velocity, **body** frame
- Translation is solved in the inertial frame and rotation in the body frame. Everything is
  referenced to the centre of mass.
- **Model selection:** `cfg.shapeModel = 'planar'` (default) or `'shape3d'`. Pass it to
  `testing/helpers/buildSeedParams.m`; `testing/helpers/seedRHS.m` then returns the matching
  right-hand side.

---

## Repository layout

```
6DOF Seed Dynamics/
├── Seed_Dynamics_ODE_Test.m       ← START HERE: a single drop, planar model
├── Seed_Dynamics_ODE_Test_3D.m    a single drop, shape3d model (+ twist/curvature, metrics)
├── README.md                      what the code is
├── RESEARCH_LOG.md                how it got here: phases, results, decisions
├── physics/                       planar model + the shared, shape-agnostic core
│   ├── seed6DOFODE.m              planar ODE right-hand side
│   ├── rigidBody6DOF.m            SHARED 6-DOF core: EOM, added-mass rate + moment, quaternions
│   ├── setupSeedShapeAndMass.m    planar builder: strips, CoM(t), inertia(t), planar switches
│   ├── validateSeedParams.m       the seed-model contract ('planar' | 'shape3d')
│   ├── translationDynamics.m      CoM linear acceleration (inertial frame, with added mass)
│   ├── rotationDynamics.m         angular acceleration (modified Euler, body frame)
│   ├── aero/
│   │   ├── computeAeroCoeffs.m    APW coefficients CT, CD, l_cp, CR ... vs angle of attack
│   │   ├── computeStripForces.m   per-strip lift + drag + rotational lift
│   │   ├── computeAngleOfAttack.m, computeStripCoP.m
│   │   ├── stripSpinDamping.m     spanwise-axis spin damping per strip (Tr) — both models
│   │   └── computeSpanForce.m, spanSpinDamping.m (Tx), normalSpinDamping.m (Ty) — planar only
│   ├── mass/                      getMassProperties.m, getAddedMass.m
│   └── helpers/                   quaternion math, per-strip velocity, Euler conversion
├── physics3d/                     the shape3d model
│   ├── setupSeedShape3D.m         builder: twist θ(s), curvature φ(s), 3D mass/inertia, switches
│   ├── seed6DOFODE3D.m            right-hand side (per-strip local frames)
│   ├── computeEdgeDrag.m          tip crossflow (edge) drag
│   ├── computeLEVForce.m          LEV vortex lift with the Rossby gate
│   └── levPlanformConstants.m     K_p, K_i, K_v for the LEV term
├── visualization/                 visualizeSeedTrajectory, visualizeSeedShape, animateSeed,
│                                  animateModeTrajectory, visualizeSeedLocalVels
├── testing/
│   ├── helpers/                   shared machinery (buildSeedParams, seedRHS, metrics, classifiers, …)
│   ├── planar/                    planar suites (+ planar/inputs/ tuned configs)
│   ├── shape3d/                   3D suites, the paper comparison, the 3D regressions
│   ├── classifier/                classifier review set builder + MODE_DEFINITIONS.md
│   ├── baselines/                 snapshot generation and model diffing
│   └── archive/                   retired scripts (Spinning_Seed_Dynamics_Test.mlx)
├── model_test_results/            baseline snapshots land here (outputs are git-ignored)
├── derivations/                   seed6DOF_physics.tex/.pdf, seed3D_physics.tex (see note)
└── Olivia Code/                   minimal_imp.m, minimal_imp_Commented.m (the 2D reference)
```

> **Derivations note.** Both `.tex` files predate the `shape3d` model and describe the
> original strip model, which is now `planar`. `seed3D_physics.tex` means "the 2D model lifted
> to 3D", not `physics3d/`. The `shape3d` physics is documented in code headers and in the
> research log until the derivation pass (roadmap).

---

## The important code

| File | Role |
|---|---|
| **[`Olivia Code/minimal_imp_Commented.m`](Olivia%20Code/minimal_imp_Commented.m)** | The **2D reference model** (Andersen–Pesavento–Wang / Pomerenk–Ristroph falling plate). The source of every sectional coefficient in both models. |
| **[`physics/rigidBody6DOF.m`](physics/rigidBody6DOF.m)** | The **shared rigid-body core**. It adds gravity, solves translation (with added mass) and modified-Euler rotation, and integrates the quaternion. It also holds both halves of Kirchhoff's added-mass pair, each behind a switch: the rate term `Ȧv` and the Munk moment `v × (A·v)`. |
| **[`physics3d/seed6DOFODE3D.m`](physics3d/seed6DOFODE3D.m)** | The **active right-hand side**. Runs the 2D APW strip aerodynamics in each strip's own frame, plus edge drag and (optionally) LEV. |
| **[`physics3d/setupSeedShape3D.m`](physics3d/setupSeedShape3D.m)** | The **shape3d builder**: twist and curvature profiles, per-strip frames, curvature-correct mass/inertia, the shape3d switch defaults. It strips the planar-only switches. |
| **[`physics/seed6DOFODE.m`](physics/seed6DOFODE.m)** | The **frozen planar right-hand side**, including the empirical whole-seed terms. |
| **[`physics/aero/computeAeroCoeffs.m`](physics/aero/computeAeroCoeffs.m)** | APW coefficients vs angle of attack, with the attached↔separated blend and three angle branches. A direct port of the 2D laws. Also holds the planar tuning constants (`C_span`, `C_span_torque`, `C_Tx`, `C_fy`, `k0_spanTorque`) and the edge-drag coefficient `C_d_edge`. |
| **[`testing/helpers/buildSeedParams.m`](testing/helpers/buildSeedParams.m)** | Builds a ready-to-run seed for either model, applying explicit switch overrides. It **warns** on a switch the chosen model does not honour instead of silently ignoring it. |

---

## Physics switches and defaults

Set a switch through `cfg` (via `buildSeedParams`) or directly on `seedParams`. The builders
own the defaults.

| switch | `planar` default | `shape3d` default | term |
|---|---|---|---|
| `enableSpanForce` | `true` | — (removed) | whole-seed spanwise force (`computeSpanForce`) |
| `enableSpanTorque` | `true` | — (removed) | moment of that force at the span CoP |
| `enableSpanCOPMigration` | `true` | — (removed) | span CoP migrates with spanwise incidence |
| `enableSpanGeomVelocity` | `false` | — (removed) | sample the span force at the geometric centre with ω×r |
| `enableSpanTorqueAttenuation` | `false` | — (removed) | reduced-frequency roll-off of the span torque |
| `enableTxDamping` | `true` | — (removed) | whole-seed roll damping `Tx` |
| `enableNormalSpinDamping` | `true` | — (removed) | whole-seed normal-axis spin damping `Ty` |
| `enableAddedMassRate` | `false` | **`true`** | `Ȧv`, translational half of Kirchhoff's pair |
| `enableAddedMassMoment` | `false` | **`true`** | `v × (A·v)` Munk moment, rotational half |
| `enableAddedMass3D` | `false` | `false` | full per-strip added-mass tensor from strip normals (meaningful for curved seeds) |
| `enableEdgeDrag` | — | `true` | tip crossflow drag, `C_d = 1.2` on `t·c̄` |
| `enableLEV` | — | `false` | LEV vortex lift (Rezgui et al. / Polhamus form) |
| `levApplicationPoint` | — | `'colocated'` | LEV acts at the APW CoP, or `'forward'` |
| `lev.rossbyDefinition` | — | `'kinematic'` | LEV gate reads `Ro = |v_ip|/(Ω·c)`, or `'geometric'` `r/c` |

Planar keeps both added-mass terms off so it stays byte-identical to every earlier run.
`shape3d` turns them on because both are first-principles and both appear in the APW reference
(eqs. 6.1–6.3). The LEV constants (`AR`, `Kp`, `Ki`, `Kv`, `RoCrit = 3`, `p = 4`,
`lambdaV = 0.5`) live on `seedParams.lev` and can be overridden with `cfg.lev`.

---

## Running the code

**Requirements:** MATLAB (R2020b or newer recommended). No toolboxes are required beyond core
MATLAB, though the Parallel Computing Toolbox speeds up the grids (`parfor`) and the Aerospace
Toolbox helps with quaternion conversion.

### Quick start: a single drop

- **Planar:** open and run **[`Seed_Dynamics_ODE_Test.m`](Seed_Dynamics_ODE_Test.m)**. It
  builds a rectangular seed with `setupSeedShapeAndMass`, integrates with `ode45`, plots the
  trajectory and Euler angles, and writes an animation.
- **shape3d:** run **[`Seed_Dynamics_ODE_Test_3D.m`](Seed_Dynamics_ODE_Test_3D.m)**. The same
  five sections, plus:
  - `twist` / `curvature` settings;
  - the shape3d switches in one block (applied through `buildSeedParams`, so a switch this
    model does not honour warns);
  - the metrics and mode label printed after the run, so one case can be checked by eye.

Scripted directly:

```matlab
addpath(genpath('physics'));  addpath('physics3d', 'visualization', 'testing/helpers');

cfg = struct('rhoFluid', 1.225, 'g', 9.81, 'shapeModel', 'shape3d');
sp  = buildSeedParams(bsp, cfg);      % bsp: seedShape, seedDensity, seedThickness,
                                      %      numStrips, tSamples, nutMass_t, nutPos_t
rhs = seedRHS(sp);                    % -> @seed6DOFODE3D here
x0  = [zeros(3,1); [1;0;0;0]; zeros(3,1); zeros(3,1)];   % level, from rest
[t, x] = ode45(@(t,x) rhs(t,x,sp), [0 5], x0);

visualizeSeedTrajectory(t, x(:,1:3).', x(:,4:7).');
```

### Test suites (`testing/`)

Each script has an editable configuration block at the top.

**shape3d**

| Script | What it does |
|---|---|
| `shape3d/runSeedTestSuite3D.m` | **The one-click 3D suite.** Edit the config block and run. Six stages, each switchable:<br>**A. Twist sweep**, released level, one animation per case.<br>**B. Curvature sweep**, bowl and one-sided, one animation per case.<br>**C. The nine 2D-parity groups** (nut mass; nut chord/span/diagonal; pitch; roll; yaw spin; asymmetry; strip convergence), using the same `seedTestCases`/`runOneSeedCase` helpers as the planar suite.<br>**D. Three moving-CoM scenarios**, animated.<br>**E. Mode grid** on our seed under the 3D model.<br>**F. Hou et al. comparison** (on by default).<br>Mode labels come from the six-part classifier (`classifySeedMode`), with the older label kept in brackets. Stages E and F report the **tracked numbers**: paper agreement, non-physical fraction, and the span tilt of revolvers. They go to `tracking_3D.txt` in the run folder, plus one row per run in `tracking_history_3D.csv` beside the run folders.<br>Physics overrides go in the `sw` struct. Stage F also writes `paper_comparison_3D.png`: ours, ours mapped to the paper's families with disagreements outlined in red, and the digitised Fig. 2a.<br>Outputs go under `model_test_results/<timestamp>_<githash>_3D_<label>/` (git-ignored). |
| `shape3d/runTwistSuite.m`, `runTwistTest.m` | Twist sweep (both-ends and single-side) with spin report; three twisted-seed animations. |
| `shape3d/runCurvatureSuite.m`, `runCurvatureTest.m` | Curvature sweep (bowl and single-side) reporting descent, cone, spin and out-of-plane CoM shift, with a flat reference row; three curved-seed animations. |
| `shape3d/exploreLEVTerm.m` | Plots the implemented LEV term against α and Rossby number, the gate, and the two application points; prints this seed's constants. |

**Paper comparison — Hou et al. (2025)**

| Script | What it does |
|---|---|
| `shape3d/runPaperAnchorCases.m` | The five configurations the paper publishes kinematics for, on their plate and release, reported against their numbers. |
| `shape3d/runPaperModeGrid.m` | Mode map over their parameter window (`x_c/a ≤ 0.416`, `y_c/b ≤ 0.208`), both classifiers per cell. Saves before reporting. |
| `shape3d/reportPaperModeGrid.m` | Reports a saved paper grid: ASCII maps, agreement with the digitised Fig. 2a, a cross-tab, and a figure. |
| `shape3d/animatePaperCases.m` | Animates the anchors and two boundary cells, to check labels by eye. |

**Classifier review** (roadmap step 1)

| File | What it does |
|---|---|
| `classifier/buildClassifierReviewSet.m` | Builds a by-eye labelling set: a coarse nut-position grid on our seed under **both** models, with run IDs shuffled so neither model nor position shows. For each run it writes:<br>• a **blind** animation (no mode colours or labels);<br>• a static summary (whole flight, attitude, cumulative rotation, rates);<br>• the raw trajectory.<br>It also writes `classifier_review.xlsx`: eye-label columns first, then the current code's labels, a blank `label_new` column, metrics, and the vocabulary. Output goes to a new timestamped folder under `Outputs/Classifier Review/`, so a labelled sheet is never overwritten. |
| `classifier/buildRelabelSheet.m` | Turns a labelled review set into a second-pass sheet with **one dropdown per classifier part**, without re-simulating:<br>• pre-filled from the first-pass labels, with only the parts each label actually states;<br>• `name_from_parts` (the classifier's rules applied to your parts) and a `MISMATCH` check;<br>• the code's answers grouped and collapsed.<br>Writes `classifier_review_v2.xlsx` next to the first pass and never modifies it. Needs Excel (driven via `actxserver`). |
| `classifier/makeModeExamples.m` | **A visual dictionary of the modes.** One idealized animation + still per mode name, plus an `INDEX.md` explaining each, written to `Outputs/Mode Examples/`. The examples are **synthetic** (prescribed textbook motion from `helpers/synthSeedMotion.m`, slowed to watchable rates), so every label has one, including modes the model never produces. Each example must classify as its own label before it is rendered. |
| `classifier/MODE_DEFINITIONS.md` | **The mode definitions**, in words and math. The six parts of the new classifier and every boundary between them, the 17 mode names in order of precedence (also the eye-labelling vocabulary), which modes are non-physical, and the two older classifiers transcribed for comparison. Kept in step with the three files below; add new modes here. |
| `helpers/classifySeedMode.m`, `defaultSeedModeThresholds.m`, `computeSeedModeMetrics.m` | **The six-part classifier.** It judges six things independently and builds the name from them: flip about the span (from attitude), revolution (from the span heading), span attitude, path, regularity, and physicality. Every boundary lives in the thresholds file, dated where it was decided on the review set. Synthetic regression: cases N1–N9 in `shape3d/testPaperMetrics.m`. |

**planar**

| Script | What it does |
|---|---|
| `planar/runSeedTestSuite.m` | Parameter sweeps (nut mass/position, initial tilt/spin, asymmetry, strip convergence); a figure per case and an overlay per group. |
| `planar/runSeedModeSuite.m` | The seven named modes at their working inputs, each classified. |
| `planar/runSeedModeGrid.m` | Chord × span nut-position phase map. **Planar only**: it calls `seed6DOFODE` directly (see roadmap). |
| `planar/pickModeGridRuns.m` | Click cells on a saved mode grid and animate them.<br>• **By default it opens the newest `model_test_results/*_3D_*/mode_grid_3D.mat`.** Set `resultsFile` to open any other grid.<br>• **Videos are written next to the grid**, in `mode_grid_picks/`, with a `picks_log.csv`.<br>• **Six-part grids get blind videos** titled with the cell's label and parts.<br>• **Model-agnostic:** it re-integrates through `seedRHS`, so planar grids work too. |
| `planar/runModeAblation.m` | The seven modes under several span-force/`Tx` configurations. |
| `planar/runSpanForceComparison.m` | Four hand-tuned cases isolating the span force, with a torque-budget diagnostic. |
| `planar/runComMovementTest.m` | Time-varying CoM (the nut slides within the body); a mode animation per scenario. |

**Regression tests**

| Script | What it checks |
|---|---|
| `shape3d/testShape3DFlatEquivalence.m` | For a flat seed, the 3D model reduces **bit-identically** to planar over the shared core (each model's extras switched off), and must *differ* at full defaults. |
| `shape3d/testCurvedMassProperties.m` | Curved-seed mass/inertia: analytic check, planar reduction, bowl sanity, toggles. |
| `shape3d/testNewPhysicsTerms.m` | The phase-4/5 terms: added-mass rate (A-checks, incl. equivalence to Kirchhoff's body-frame form), Munk moment (M-checks, incl. the 2D limit vs APW eq. 6.3), edge drag, and LEV (constants, both Rossby gates, application points, RHS wiring, `cfg.lev` overrides). |
| `shape3d/testPaperMetrics.m` | The paper metrics and classifier on **synthetic** motion with known answers (8 cases), including the tilted-revolution trap. |

**Baselines**

| Script | What it does |
|---|---|
| `baselines/generateModelBaseline.m` | Timestamped, git-stamped snapshot into `model_test_results/`. It runs these stages:<br>• the coarse mode grid, CoM scenarios and test suite, under **both** models on the same flat seed (`mode_grid/` vs `mode_grid_3d/`, etc.);<br>• a `shape3d/` twist + curvature stage;<br>• the `paper/` comparison (anchors + grid), regenerated every time.<br>It records both models' switch sets. |
| `baselines/compareBaselineModels.m` | Diffs a snapshot's planar and shape3d stages: mode census, migration table, metric shifts, a changed-cell map. Optionally regression-checks planar against an older snapshot. |

---

## Reproducing the flight modes (planar model)

The configuration below reproduces the seven reference modes in the **planar** model and is
its builder default ("FULL-minus-geomVelocity"). It is recorded in
[`testing/planar/inputs/Working Dynamics Inputs 8-2-26.txt`](testing/planar/inputs/Working%20Dynamics%20Inputs%208-2-26.txt).
The shape3d model does not use these terms; see the switch table above.

**Base seed:**

- nut mass `75e-6` kg at the body centre;
- body density `65` kg/m³, thickness `0.002` m;
- span `0.050` m, chord `0.015` m;
- air (`rhoFluid = 1.225`, buoyancy ignored);
- released from rest for 10 s unless an initial attitude is given.

**Planar switches:** the planar column of the switch table. **Tuned aero constants**
(`computeAeroCoeffs` defaults): `C_span = 0.2`, `C_span_torque = 0.7`, `C_Tx = 1.0`,
`C_fy = 1`, `k0_spanTorque = 0.2`. The attenuation is not needed with this tuning.

**Mode-eliciting inputs** (`c = chordLength`, `S = spanLength`; `theta0 = [0 0 pi/6]` is a
π/6 tilt about body `z`):

| Mode | Nut position `[x; y; z]` | Initial condition |
|---|---|---|
| Spanwise-axis fluttering | `[0; 0; 0]` (center) | `theta0 = [0 0 pi/6]` |
| Gliding | `[0.5*c; 0; 0]` | from rest |
| Diving | `[1.5*c; 0; 0]` | from rest |
| Fluttering + spiral | `[0; 0; 0.01*S]` | `theta0 = [0 0 pi/6]` |
| Fluttering + tight spiral | `[0; 0; S]` | `theta0 = [0 0 pi/6]` |
| Autorotation | `[c; 0; 1.2*S]` | from rest |
| Parachute | `[0; -c; 0]` | from rest |

> Illustrative examples, not mode-onset boundaries. Several of these nut positions lie
> outside the planform (e.g. autorotation at 1.2 S), which no real samara has.

---

## The model in brief

**Both models**

- **Per strip:**
  - translational lift and drag at the migrating APW centre of pressure;
  - rotational lift at mid-chord;
  - spanwise-axis spin damping (`Tr`).
  
  All from the 2D quasi-steady coefficient laws.
- **Rigid body:** gravity, added mass (flat-plate form from per-strip 2D added mass),
  gyroscopic and time-varying-inertia torques.
- **Omitted:** buoyancy (negligible in air) and the rotational added-inertia rate `Ȧ_rot·ω`.

**planar only.** The empirical whole-seed terms: the span force and its torque at a migrating
span CoP, `Tx`, `Ty`.

**shape3d only**

- twist and curvature, with each strip's aero in its own frame;
- curvature-correct mass/inertia;
- Kirchhoff's added-mass pair (rate + Munk moment), on by default;
- tip edge drag, on by default;
- LEV vortex lift with a kinematic Rossby gate, off by default.

---

## References

**Reference model and falling plates**

1. Pomerenk, O., & Ristroph, L. (2024). *Aerodynamic equilibria and flight stability of
   plates at intermediate Reynolds numbers.* J. Fluid Mech. [arXiv:2408.08864](https://arxiv.org/abs/2408.08864).
   The steady equilibria (gliding, diving) and stability of falling plates; the basis of the
   2D reference code.
2. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Unsteady aerodynamics of fluttering and
   tumbling plates.* J. Fluid Mech. 541, 65–90. The 2D quasi-steady model this code extends:
   the sectional coefficients, and eqs. (6.1)–(6.4) for added mass and the `(m₁₁ − m₂₂)`
   moment.
3. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Analysis of transitions between
   fluttering, tumbling and steady descent of falling cards.* J. Fluid Mech. 541, 91–104.
   The flutter/tumble transition the 2D model reproduces.

**Samaras: experiments and target data**

4. Hou, Z.-B., Zhang, J.-D., Li, Y.-D., Jia, Y.-X., & Huang, W.-X. (2025). *Aerodynamic
   significance of mass distribution on diverse samara descent behaviors.* Commun. Eng. 4,
   129. [doi:10.1038/s44172-025-00465-8](https://doi.org/10.1038/s44172-025-00465-8).
   Experimental mode map (AR / ST / CH / FA) over CoM position. Its Fig. 2a and five anchor
   cases are the phase-5 validation target.
5. Lentink, D., Dickson, W. B., van Leeuwen, J. L., & Dickinson, M. H. (2009). *Leading-edge
   vortices elevate lift of autorotating plant seeds.* Science 324, 1438–1440. Shows the LEV
   that strip theory omits, and its lift levels (~1.8–2.0), which the LEV term matches.
6. *Mechanism of autorotation flight of maple samaras (Acer palmatum).* (2014). Exp. Fluids.
   [doi:10.1007/s00348-014-1718-4](https://doi.org/10.1007/s00348-014-1718-4). Measured
   autorotation kinematics (descent speed, spin rate, coning angle).
7. Norberg, R. Å. (1973). *Autorotation, self-stability, and structure of single-winged
   fruits and seeds (samaras) with comparative remarks on animal flight.* Biol. Rev. 48,
   561–596. Cited for the chordwise CoM position (27–35% chord) giving pitch stability,
   against the model's chordwise sweep (log, phase 4).

**Leading-edge vortex and vortex lift**

8. Lentink, D., & Dickinson, M. H. (2009). *Rotational accelerations stabilize leading edge
   vortices on revolving fly wings.* J. Exp. Biol. 212, 2705–2719. The Rossby-number
   dependence of LEV stability behind the LEV gate (no critical value is published).
9. Rossby, C.-G. (1936). *Dynamics of steady ocean currents in the light of experimental
   fluid mechanics.* Papers in Physical Oceanography and Meteorology 5(1). The general
   Rossby ratio `U/(Ω·L)`, used for the kinematic gate definition.
10. Rezgui, D., Arroyo, J. C., & Theunissen, R. (2020). The Aeronautical Journal 124(1278),
    1236–1261. [doi:10.1017/aer.2020.25](https://doi.org/10.1017/aer.2020.25). Sectional
    adaptation of the Polhamus analogy to a samara blade; the source of the LEV term's form
    (eqs. 3–5).
11. Polhamus, E. C. (1966). *A concept of the vortex lift of sharp-edge delta wings based on
    a leading-edge-suction analogy.* NASA TN D-3767. The original suction analogy:
    potential lift + vortex lift.
12. Polhamus, E. C. (1971). *Predictions of vortex-lift characteristics by a leading-edge
    suction analogy.* J. Aircraft 8(4), 193–199. The generalised form of the analogy used by
    [10].
13. Snyder, M. H., & Lamar, J. E. (1972). *Application of the leading-edge-suction analogy to
    prediction of longitudinal load distribution and pitching moments for sharp-edged delta
    wings.* NASA TN D-6994. Vortex-lift loading sits near the potential-flow CoP; the basis
    for the `'colocated'` LEV application point.
14. Lamar, J. E. (1974). *Extension of leading-edge-suction analogy to wings with separated
    flow around the side edges at subsonic speeds.* NASA TR R-428. Side-edge (tip) vortex
    lift; planned low-AR term (roadmap 3d).

**Finite-wing theory**

15. Prandtl, L. (1918). *Tragflügeltheorie.* Nachr. Ges. Wiss. Göttingen. Lifting-line
    theory and the elliptic induced-drag factor used for `K_i`.
16. Helmbold, H. B. (1942). *Der unverwundene Ellipsenflügel als tragende Wirbelfläche.*
    Jahrbuch der Deutschen Luftfahrtforschung. Low-aspect-ratio lift slope, used for `K_p`
    and for the planned finite-AR lift correction.
17. Jones, R. T. (1946). *Properties of low-aspect-ratio pointed wings at speeds below and
    above the speed of sound.* NACA Report 835. Slender-wing lift. It is the yardstick the
    planar span patch was checked against (log, phase 6) and the planned term for
    near-spanwise flow.
18. Glauert, H. (1935). *Airplane propellers.* In W. F. Durand (ed.), *Aerodynamic Theory*,
    Vol. IV. Momentum theory, the windmill-brake state and the Prandtl tip-loss factor, used
    to rule out induced inflow and tip loss (log, phase 4).

**Rigid bodies in fluid**

19. Lamb, H. (1932). *Hydrodynamics* (6th ed.), ch. VI. Cambridge University Press.
    Kirchhoff's equations for a body in fluid: the added-mass rate and the Munk moment.

**Canonical static plate data (planned calibration, roadmap step 2)**

20. Taira, K., & Colonius, T. (2009). *Three-dimensional flows around low-aspect-ratio
    flat-plate wings at low Reynolds numbers.* J. Fluid Mech. 623, 187–207. Low-AR plate
    forces and wake structure at Re 300–500, post-stall.
21. Ortiz, X., Rival, D., & Wood, D. (2015). *Forces and moments on flat plates of small
    aspect ratio with application to PV wind loads and small wind turbine blades.* Energies
    8(4), 2438–2453. [doi:10.3390/en8042438](https://doi.org/10.3390/en8042438). Tunnel
    forces and moments for AR 0.4–9 plates, Re 6×10⁴–2×10⁵.
22. Shields, M., & Mohseni, K. (2012). *Effects of sideslip on the aerodynamics of
    low-aspect-ratio low-Reynolds-number wings.* AIAA J. 50(1), 85–99.
    [doi:10.2514/1.J051151](https://doi.org/10.2514/1.J051151). Sideslip effects on low-AR
    wings, the physics behind the planned sideslip roll moment.

---

## Roadmap

### Next — decided 2026-09-26 (rationale: [log, phase 6](RESEARCH_LOG.md#phase-6--3d-test-suite-and-the-rollslide-diagnosis-2026-09-24--09-26))

**Rule for any new term.** Its form comes from theory, and its magnitude from canonical
*static* published data. Free-flight data (Hou et al., our own trials) is validation only,
never tuning.

- [ ] **1. Fix the classifier.** Every later step is judged through it.
  - Add explicit **non-physical labels** (end-on / span-down, edge-slide). They count as
    disagreements, and their grid fraction is a tracked metric that should go to zero.
  - Build a **hand-labelled reference set** of saved trajectories labelled by eye. Tune
    thresholds on half, check on the other half, then freeze it as a regression test on the
    *saved* trajectories.
  - [x] Tooling: `classifier/buildClassifierReviewSet.m` and `classifier/MODE_DEFINITIONS.md`.
  - [ ] Label the set (first pass done 2026-10-03).
  - [ ] Paper-metric sign error fixed: steady revolvers no longer read as spiral tumbling.
  - [ ] Six-part classifier written (`classifySeedMode`), synthetic regression passing.
  - [x] Visual dictionary of the modes (`classifier/makeModeExamples.m` → `Outputs/Mode Examples`).
  - [ ] Re-run the paper comparison: every paper-classifier result before 2026-10-03 has the
        sign error.
- [ ] **2. Virtual wind tunnel.** Evaluate the aero model on a plate held fixed at any
  orientation, with imposed rotation rates and no integration: lift, drag, chordwise and
  spanwise CoP, rolling moment vs sideslip, the end-on moment, roll damping.
  - Compare against refs [20]–[22].
  - Run it on the planar model too, to see what each patch contributes.
  - This subsumes the open phase-5 question: the static pitch moment vs chordwise CoM, i.e.
    why `y_c/b` gives so little pitch stiffness.
- [ ] **3. Low-AR physics, one piece at a time**, each accepted only if it improves step 2:
  - [ ] **a. Sideslip roll moment.** Redistribute the existing strip normal loads toward the
    windward tip as a function of local spanwise incidence, with zero net force added. The
    principled successor to the planar span torque.
  - [ ] **b. Slender-wing lift for near-spanwise flow.** Jones [17] at small angles, continued
    with the Polhamus/Lamar potential + vortex form. Removes the span-down mode.
  - [ ] **c. Finite-AR attached lift** (Helmbold [16]). APW's `CL1 = 5.2` vs 3.4–3.6/rad at
    AR 3–3.3. *(This is the planform correction, distinct from the rotor tip-loss factor
    dropped in phase 4.)*
  - [ ] **d. Side-edge vortex lift** (Lamar [14]), later.
- [ ] **4. Validation protocol for our seed trials.**
  - Measure shape, mass, CoM, mass distribution and release conditions.
  - Record model predictions *before* viewing the flights.
  - Repeat drops per seed to measure scatter.
  - Monte-Carlo the model over the CoM and release uncertainty, so both sides give mode
    probabilities rather than single labels. Apply the same treatment to Hou et al.

### Other open items

- [ ] **Second-pass relabelling (deferred 2026-10-04).** The main classifier issues are fixed;
  scoring against eye labels can wait. Everything is ready:
  - `classifier_review_v2.xlsx` in the review set: one dropdown per part, pre-filled, with a
    `MISMATCH` check (`classifier/buildRelabelSheet.m`).
  - The steps: relabel the flagged runs (MODE_DEFINITIONS section 3), then score
    `label_new` on a held-out half.

- [ ] `planar/runSeedModeGrid.m` hardcodes `seed6DOFODE` instead of `seedRHS`, so it cannot run
  the 3D model. The 3D suite's stage E is the 3D grid in the meantime.
- [ ] Two-point-mass support (`bsp.extraMasses`, shape3d builder only), to carry the paper
  plate's second weight rather than one equivalent nut. That nut is ~1.6× low in long-axis
  inertia at the top of the `y_c/b` range.
- [ ] Nut form drag. The nut is a drag-free point mass, though in autorotation it sweeps
  off-axis.
- [ ] Re-measure our seed's descent against `sqrt(σg/ρ)` under the current shape3d defaults.
  The old "autorotation descends too fast" bug predates the Kirchhoff defaults.
- [ ] 2D validation suite vs `minimal_imp`: RHS match at matched states, plus trajectory
  comparison.
- [ ] Standalone intermediates-recovery helper (rebuild per-strip quantities over a saved
  trajectory; `ode45` discards the RHS's second output).
- [ ] Validate twist and curvature against real seeds; optionally accept a general YZ polyline
  profile instead of the parametric twist/curvature functions.
- [ ] Tapered / arbitrary planform support beyond the rectangular seed.
- [ ] Momentum-rigorous internal mass movement (the moving-nut test prescribes the CoM
  kinematically; there is no reaction from the sliding mass).
- [ ] **Derivation pass** (once the model is finalized):
  - Fix `seed6DOF_physics.tex`. It says `Tx` is off by default (it is on in planar) and that
    `Ty` is always applied (it has a switch). Refresh its PDF.
  - Rename `seed3D_physics.tex` so it cannot be confused with `physics3d/`.
  - Write a `shape3d` derivation. First candidates, since they are settled: twist/curvature
    frames, curved mass properties, edge drag, and Kirchhoff's added-mass pair.
- [ ] Housekeeping:
  - Stale code comments: `seed6DOFODE3D.m`'s header and `runPaperAnchorCases.m`'s header still
    describe the added-mass rate as off.
  - MATLAB autosave files (`*.asv`) are tracked; add them to `.gitignore`.

### Completed

Details and measurements for each are in the [research log](RESEARCH_LOG.md).

- [x] Flat-plate 6-DOF strip model reproducing the seven reference modes (planar; phase 0).
- [x] Automatic flight-mode classification and the chord × span phase map (phase 0).
  Calibrated on the reference inputs; now under review (roadmap step 1).
- [x] Outputs: trajectory + Euler plots, seed animations, sweep figures, trajectory metrics,
  torque-budget diagnostic, moving-CoM test, mode animations.
- [x] Compiled PDF of `seed6DOF_physics.tex`. It predates later model changes; see the
  derivation pass.
- [x] Non-planar seed shape: the `shape3d` model with twist, curvature and curved mass
  properties (phase 1).
- [x] Phase 1: curvature benchmark.
- [x] Phase 2: invented terms cut from `shape3d`. Resolved there are:
  - the `Tx` double count;
  - the phantom span torque;
  - the dimensionally wrong `Ty`.
  
  Planar keeps all three by design.
- [x] Phase 3: benchmark gate.
- [x] Phase 4: real physics back.
  - Edge crossflow drag implemented (replaces the span force).
  - Added-mass rate `-Ȧv` implemented.
  - LEV vortex lift implemented (off by default).
  - Nut form drag remains open (above).
- [x] Phase 5, Tier 0: paper-parameter mode grid, paper-mode classifier, anchor cases, and the
  synthetic-motion regression for the metrics.
- [x] Phase 5, Tier 1:
  - Added-mass rate on by default in shape3d.
  - Munk moment implemented and on by default.
  - LEV screened: it stays off.
- [x] One-click 3D test suite and single-drop 3D harness (phase 6).

### Not doing

These were considered and dropped or superseded. Reasons are in the
[log's register](RESEARCH_LOG.md#decided-against-or-superseded).

- Tuning the span torque.
- The reduced-frequency attenuation and its unsteady-CoP-lag replacement.
- A crossflow-drag-only span-force variant.
- Induced inflow and Prandtl rotor tip loss.
- LEV as the collapse fix.
- Moving edge drag to the windward tip.
- Pruning the planar constants.
- Fixing the `(π−|α|)` branch quirks.
- The "combined" two-plane re-discretization (set aside; its problem is revisited by step 3).
