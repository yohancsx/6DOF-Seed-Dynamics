# 6-DOF Seed Dynamics
Create by Yohan Sequeira, with great help from Olivia Pomerenk

A MATLAB simulation of the free flight of an autorotating winged seed (a **samara**) —
modelled as a flat plate with a discrete "nut" mass — using **quasi-steady strip-theory
aerodynamics**. It is a full 6-degree-of-freedom (13-state) extension of the classic 2D
Andersen–Pesavento–Wang/Pomerenk-Ristroph falling-plate model into three dimensions.

The goal is to reproduce, and let you explore, the diverse descent behaviours real seeds
exhibit — gliding, diving, fluttering, tumbling/spiralling, and autorotation — and to see
how they emerge from **where the seed's mass sits** relative to its wing.

---

## What it does

The wing is sliced into spanwise **strips**; each strip is treated as a local 2D
cross-section using the exact 2D coefficient laws, and the per-strip forces and torques are
summed and integrated with `ode45`. Because it's a faithful 3D lift of the 2D model, a
single strip confined to the plane reproduces the validated 2D falling-plate dynamics, while
the full 3D model adds quaternion attitude, gyroscopic and time-varying-inertia torques, and
optional whole-seed spanwise effects.

Move the nut around and you change the centre of mass, which is the primary control on the
flight mode — the same lever real samaras use.

---

## Conventions (worth knowing before reading the code)

- **Body axes:** `x` = chordwise, `y` = plate-normal, `z` = spanwise. The plate lies in the
  body `x`–`z` plane.
- **World frame:** Y-up. Gravity is `[0; -g; 0]`.
- **State vector (13×1):** `x = [ r(3) ; q(4) ; v(3) ; omega(3) ]`
  - `r` — CoM position, **inertial** frame
  - `q` — orientation quaternion `[q0;q1;q2;q3]`, scalar-first, **body → world**
  - `v` — CoM velocity, **inertial** frame
  - `omega` — angular velocity, **body** frame
- Translation is solved in the inertial frame; rotation in the body frame; everything is
  referenced to the centre of mass.

The full physics — every force, torque, aerodynamic-coefficient branch, and added-mass term
— is derived in **[`derivations/seed6DOF_physics.tex`](derivations/seed6DOF_physics.tex)**.
Read that for the math; this README is orientation.

---

## Repository layout

```
6DOF Seed Dynamics/
├── Seed_Dynamics_ODE_Test.m       ← START HERE: script that runs a simple drop
├── physics/                       planar model + the shared, shape-agnostic core
│   ├── seed6DOFODE.m              planar ODE right-hand side (the integrator function)
│   ├── rigidBody6DOF.m            SHARED 6-DOF core: EOM + quaternion kinematics
│   ├── setupSeedShapeAndMass.m    builds the planar seed: strips, CoM(t), inertia(t)
│   ├── validateSeedParams.m       the seed-model contract ('planar' | 'shape3d')
│   ├── translationDynamics.m      CoM linear acceleration (inertial frame)
│   ├── rotationDynamics.m         angular acceleration (modified Euler, body frame)
│   ├── aero/                      aerodynamics
│   │   ├── computeAeroCoeffs.m    CT, CD, l_cp, CR ... vs angle of attack
│   │   ├── computeStripForces.m   per-strip lift + drag + rotational lift
│   │   ├── computeAngleOfAttack.m, computeStripCoP.m
│   │   ├── stripSpinDamping.m     spanwise-axis spin damping (Tr)
│   │   ├── spanSpinDamping.m, normalSpinDamping.m, computeSpanForce.m   (whole-seed)
│   │   └── ...
│   ├── mass/                      getMassProperties.m, getAddedMass.m
│   └── helpers/                   quaternion math, per-strip velocity, Euler conversion
├── physics3d/                     NON-PLANAR model (twist / curvature) — see the roadmap
│   ├── setupSeedShape3D.m         3D builder: twist θ(s), curvature φ(s), 3D mass/inertia
│   └── seed6DOFODE3D.m            3D right-hand side (per-strip local frames)
├── visualization/                 visualizeSeedTrajectory.m, animateModeTrajectory.m, ...
├── testing/
│   ├── helpers/                   shared machinery (buildSeedParams, seedRHS, classifier, …)
│   ├── planar/                    planar suites (+ planar/inputs/ reference configs)
│   ├── shape3d/                   twist + curvature suites, the 3D regressions
│   ├── baselines/                 generateModelBaseline.m (snapshot tooling)
│   └── archive/                   retired scripts
├── model_test_results/            baseline snapshots (outputs are git-ignored)
├── derivations/                   seed6DOF_physics.tex   (the physics writeup)
└── Olivia Code/                   minimal_imp_Commented.m   (the 2D reference model)
```

---

## The important code

| File | Role |
|---|---|
| **[`Olivia Code/minimal_imp_Commented.m`](Olivia%20Code/minimal_imp_Commented.m)** | The **2D reference model** (Andersen–Pesavento–Wang//Pomerenk-Ristroph falling plate). Ground truth for all aerodynamic physics — every coefficient and force term in the 3D code traces back to it. |
| **[`physics/seed6DOFODE.m`](physics/seed6DOFODE.m)** | The **ODE right-hand side** passed to `ode45`. Orchestrates mass properties, per-strip aero, whole-seed terms, and the equations of motion; returns the 13-state derivative. This is the heart of the simulation. |
| **[`physics/aero/computeAeroCoeffs.m`](physics/aero/computeAeroCoeffs.m)** | The **aerodynamic coefficients** (`CT`, `CD`, centre-of-pressure fraction, rotational-lift and spin-damping constants) as functions of angle of attack, with the attached↔separated blend and the three angle-of-attack branches. A direct port of the 2D coefficient laws. |
| **[`physics/setupSeedShapeAndMass.m`](physics/setupSeedShapeAndMass.m)** | Turns a 2D wing polyshape + a nut mass into a full `seedParams` struct: strip geometry, and the (optionally time-varying) CoM and inertia tensor. |
| **[`physics/rigidBody6DOF.m`](physics/rigidBody6DOF.m)** | The **shape-agnostic rigid-body core**: given mass properties and the net body-frame force/torque, it adds gravity, solves the translational (with added mass) and modified-Euler rotational dynamics, and integrates the quaternion. Both the planar and the 3D models call it, so the dynamics that must never differ live in one place. |
| **[`physics3d/setupSeedShape3D.m`](physics3d/setupSeedShape3D.m)** | The **non-planar builder**: spanwise twist and curvature, per-strip local frames, and curvature-correct mass/inertia. Paired with [`physics3d/seed6DOFODE3D.m`](physics3d/seed6DOFODE3D.m), which runs the same 2D strip aero in each strip's own frame. Selected via `cfg.shapeModel = 'shape3d'`. |

> **Note on the spanwise-flow additions.** `computeSpanForce`, `spanSpinDamping` (`Tx`),
> and `normalSpinDamping` (`Ty`) are experimental whole-seed extensions gated behind
> `seedParams.enable*` switches. The validated, 2D-reducible baseline runs with the
> spanwise force off; see the derivation for what each switch does.

---

## Running the code

**Requirements:** MATLAB (R2020b or newer recommended). No toolboxes required beyond core
MATLAB. But it's nice if you have the aerospace toolbox for quaternion conversion.

### Quick start — a single drop

Open and run **[`Seed_Dynamics_ODE_Test.m`](Seed_Dynamics_ODE_Test.m)**. It:

1. adds `physics/` and `visualization/` to the path,
2. builds a simple rectangular seed via `setupSeedShapeAndMass`,
3. sets an initial state and integrates with `ode45`, and
4. plots the 3D trajectory and Euler-angle history.

Edit the seed geometry, nut position, and initial conditions at the top to explore.

If you'd rather script it directly:

```matlab
addpath(genpath('physics'));
addpath('visualization');

seedParams = setupSeedShapeAndMass(seedParamsIn);   % build the seed
seedParams.rhoFluid = 1.225;   seedParams.g = 9.81; % air, gravity

x0 = [zeros(3,1); [1;0;0;0]; [0;-0.3;0]; [0;0;20]]; % 13-state: level, small fall + spin
[t, x] = ode45(@(t,x) seed6DOFODE(t,x,seedParams), [0 2], x0);

visualizeSeedTrajectory(t, x(:,1:3).', x(:,4:7).');
```

### Test suites (`testing/`)

`testing/` is split by model: `planar/` for the flat-plate suites, `shape3d/` for the
non-planar ones, `helpers/` for machinery both share, and `baselines/` for snapshots.

| Script | What it does |
|---|---|
| `planar/runSeedTestSuite.m` | Parameter sweeps (nut mass/position, initial tilt/spin, asymmetry, strip-count convergence); saves a trajectory figure per case and a per-group overlay. |
| `planar/runSeedModeSuite.m` | Elicits and auto-classifies the biological flight modes (glide, dive, spiral, autorotation, parachute). |
| `planar/runSeedModeGrid.m` | Sweeps a chord×span grid of nut positions and renders the flight-mode phase map. `planar/pickModeGridRuns.m` lets you click cells on that map and animate them. |
| `planar/runSpanForceComparison.m` | Four hand-tuned cases isolating the spanwise-force physics, with a torque-budget diagnostic. |
| `planar/runComMovementTest.m` | Drives a **time-varying** CoM (the nut slides within the body) to test mode transitions; renders a mode animation per scenario. |
| `shape3d/runTwistSuite.m` | Sweeps spanwise twist (both-ends and single-side) and reports the resulting spin. |
| `shape3d/runTwistTest.m` | Renders animations for three twisted seeds (slight / large / asymmetric twist). |
| `shape3d/runCurvatureSuite.m` | Sweeps spanwise curvature (symmetric bowl and single-side) and reports descent, cone, spin and the out-of-plane CoM shift. Includes a **flat reference row**: a zero-valued curvature profile still takes the curved branch, which switches the span force off, so the 0° row is not the flat baseline. |
| `shape3d/runCurvatureTest.m` | Renders animations for three curved seeds (slight / large / asymmetric curvature). |
| `baselines/generateModelBaseline.m` | Timestamped, git-stamped snapshot into `model_test_results/`. Runs the coarse mode grid, CoM-movement scenarios and sweep test suite under **both models on the same flat seed** (`mode_grid/` vs `mode_grid_3d/`, etc.), plus a `shape3d/` stage sweeping twist and curvature. Records both models' switch sets. |
| `baselines/compareBaselineModels.m` | Diffs a snapshot's planar and shape3d stages: mode census, migration table (which mode became which), metric shifts on unchanged cells, a three-panel map with a changed-cell overlay, and the CoM/test-suite outcomes. Optionally regression-checks the planar stages against an older snapshot. **Re-run after every phase-4 change.** |
| `shape3d/testShape3DFlatEquivalence.m` | Regression: the 3D model must reduce **bit-identically** to the planar one for a flat seed. |
| `shape3d/testCurvedMassProperties.m` | Regression: curved-seed mass/inertia (analytic check, planar reduction, bowl sanity, toggles). |
| `baselines/generateModelBaseline.m` | Writes a timestamped, git-stamped snapshot (coarse mode grid + CoM suite + test suite + a config manifest) to `model_test_results/`. |

Each script has an editable configuration block at the top. The suites add the paths they
need; a bare session wants `physics/`, `physics3d/`, `visualization/` and `testing/helpers/`.

---

## Reproducing the flight modes

The configuration below reproduces the full set of biological descent modes. **It is now
the code default** ("FULL-minus-geomVelocity"), set in `setupSeedShapeAndMass` /
`buildSeedParams` / `seed6DOFODE`, so a seed built the normal way already uses it — you do
not need to set these by hand. Recorded in
[`testing/planar/inputs/Working Dynamics Inputs 8-2-26.txt`](testing/planar/inputs/Working%20Dynamics%20Inputs%208-2-26.txt).

**Base seed** (all runs): nut mass `75e-6` kg at the body center, body density `65` kg/m³,
span `0.050` m, chord `0.015` m; air (`rhoFluid = 1.225`, buoyancy ignored); released from
rest for 10 s unless an initial attitude is given.

**Physics switches** (`seedParams.*`, the defaults):

| switch | value |
|---|---|
| `enableSpanForce` | `true` |
| `enableSpanTorque` | `true` |
| `enableSpanGeomVelocity` | `false` |
| `enableSpanCOPMigration` | `true` |
| `enableSpanTorqueAttenuation` | `false` |
| `enableTxDamping` | `true` |

**Tuned aero** (`seedParams.aero.*`; every other coefficient stays at its default):
`C_span = 0.2`, `C_span_torque = 0.7` (with `C_fy = 1`, `C_Tx = 1.0`, `k0_spanTorque = 0.2`).
With this tuning the reduced-frequency torque attenuation is *not* needed — stability comes
from the weak span force plus roll damping (`Tx`).

**Mode-eliciting inputs** (`c = chordLength`, `S = spanLength`; `theta0 = [0 0 pi/6]` is an
initial tilt of π/6 about the body-`z` / spanwise axis):

| Mode | Nut position `[x; y; z]` | Initial condition |
|---|---|---|
| Spanwise-axis fluttering | `[0; 0; 0]` (center) | `theta0 = [0 0 pi/6]` |
| Gliding | `[0.5*c; 0; 0]` | from rest |
| Diving | `[1.5*c; 0; 0]` | from rest |
| Fluttering + spiral | `[0; 0; 0.01*S]` | `theta0 = [0 0 pi/6]` |
| Fluttering + tight spiral | `[0; 0; S]` | `theta0 = [0 0 pi/6]` |
| Autorotation | `[c; 0; 1.2*S]` | from rest |
| Parachute | `[0; -c; 0]` | from rest |

> These are basic illustrative examples — not the exact mode-onset boundaries, and not
> necessarily the cleanest instance of each mode.

---

## The model in brief

- **Per strip:** translational lift + drag (at the migrating centre of pressure) and
  rotational lift (at mid-chord), from the 2D quasi-steady coefficient laws.
- **Whole seed:** gravity, added mass, gyroscopic and time-varying-inertia damping torques, and the
  optional spanwise-flow force/torque and spin-damping terms.
- **Deliberately omitted:** buoyancy (negligible in air) and added-mass CoM-offset coupling.

Full derivation, including every coefficient branch and added-mass term, is in
[`derivations/seed6DOF_physics.tex`](derivations/seed6DOF_physics.tex).

---

## References

The model and its validation targets draw on (not an exhaustive list):

1. Pomerenk, O., & Ristroph, L. (2024). *Aerodynamic equilibria and flight stability of
   plates at intermediate Reynolds numbers.* Journal of Fluid Mechanics.
   [arXiv:2408.08864](https://arxiv.org/abs/2408.08864) — steady equilibria (gliding,
   diving) and stability of the falling-plate modes.
2. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Unsteady aerodynamics of fluttering
   and tumbling plates.* Journal of Fluid Mechanics, 541, 65–90. — the 2D quasi-steady
   falling-plate model this code extends.
3. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Analysis of transitions between
   fluttering, tumbling and steady descent of falling cards.* Journal of Fluid Mechanics,
   541, 91–104.
4. *Aerodynamic significance of mass distribution on diverse samara descent behaviors.*
   (2025). Communications Engineering.
   [nature.com/articles/s44172-025-00465-8](https://www.nature.com/articles/s44172-025-00465-8)
   — how CoM position maps to samara flight modes (autorotation, spiral tumbling, chaotic,
   falling).
5. Lentink, D., Dickson, W. B., van Leeuwen, J. L., & Dickinson, M. H. (2009).
   *Leading-edge vortices elevate lift of autorotating plant seeds.* Science, 324, 1438–1440.
   — the leading-edge vortex that a strip-theory model does not capture.
6. *Mechanism of autorotation flight of maple samaras (Acer palmatum).* (2014). Experiments
   in Fluids. [doi:10.1007/s00348-014-1718-4](https://doi.org/10.1007/s00348-014-1718-4)
   — measured autorotation kinematics (descent speed, spin rate, coning angle).

---

## TODO / roadmap

Working checklist of what's outstanding. Check items off (`- [x]`) as they land.

### Model cleanup — cut the invented 3D terms (active arc)

A physics audit (Sep 2026) found that every invented term in the model exists to compensate
for one structural limitation: a flat, single-plane strip decomposition cannot produce
out-of-plane force or moment. **Twist and curvature now supply that from real geometry**, so
most of those terms are scaffolding rather than physics. The measurements behind this are
recorded as individual items under *Bugs* below.

Goal: get back to **Andersen–Pesavento–Wang + strips + first principles + added mass**, then
measure what survives before adding anything new. The phases run in order; phase 3 gates phase 4.

- [ ] **1 — Curvature tests in the benchmark suite.** *(built and run — confirm before ticking.)*
  `testing/shape3d/runCurvatureSuite.m` (symmetric bowl + single-side, swept 0–35° tip dihedral,
  with a separate **flat reference row**) and `runCurvatureTest.m` (three animations). Baseline
  captured before the cuts. Result: a **symmetric bowl produces no spin at all** (0.00 rad/s
  throughout) but collapses the cone 51°→0.1° and turns the seed into a parachute; **asymmetric**
  curvature does spin it (up to ~50 rad/s). The out-of-plane CoM shift is curvature's signature
  (+1.40 mm at 35°, matching `testCurvedMassProperties`'s hand calculation exactly).
  Note the flat-vs-0° caveat: a zero-valued curvature profile still takes the curved branch and
  switches the span force off, so the 0° row is not the flat baseline. Measured, that config
  difference is nil here (1.66 vs 1.66 m/s descent, 51.0° vs 50.9° cone) — itself consistent with
  the span force being 1% of weight.
- [ ] **2 — Cut the invented terms, in `physics3d/` only.** Remove `Tx`, `Ty`, the span torque,
  the span-CoP migration and the reduced-frequency attenuation from `seed6DOFODE3D.m`, along
  with the constants `C_span_torque` / `k0_spanTorque` / `C_Tx` / `C_fy`. The planar RHS stays
  frozen and bit-identical, and the shared `computeAeroCoeffs` is untouched — the 3D RHS simply
  stops reading the dead constants. This is what the dual codebase was built for.
  *Note:* the span **force** (`computeSpanForce` / `C_span`) is deliberately **not** cut here;
  only its torque apparatus is. Its moment becomes the honest `cross(r, F)` at the real
  application point (measured `1.07e-07` N·m — negligible, but correct).

  *(done — confirm before ticking.)* The switch set now differs by model: a `shape3d`
  seedParams carries only `enableSpanForce`, `enableSpanGeomVelocity`, `enableAddedMass3D`,
  because `setupSeedShape3D` strips the rest and `buildSeedParams` **warns** on a stale
  `cfg.enable*` or `cfg.aero.*` override instead of applying it silently. One prerequisite was
  added to the planar model: `enableNormalSpinDamping` (default `true`, so planar behaviour is
  unchanged) — `Ty` was previously unconditional, leaving no way to configure planar down to the
  bare core. Verified: all seven planar modes still classify correctly, and
  `testShape3DFlatEquivalence` was **re-scoped rather than deleted** — it now asserts the two
  models match **bit-for-bit over the shared core** (`0.000e+00`, fixed and moving nut) *and*
  that they still differ at full defaults, so a cut term silently creeping back would fail it.

  **What the cuts changed.** Curvature: nothing measurable — a bowl's ω stays near zero, so
  `Tx ∝ ω²` was already inert there, and `enableSpanForce` was already off for curved seeds.
  Twist: substantially. Anti-symmetric family, pre-cut → post-cut:

  | tip twist | 0° | 5° | 10° | 15° | 20° | 25° | 30° |
  |---|---|---|---|---|---|---|---|
  | spin, pre (rad/s) | 0.0 | 40.1 | 56.2 | 61.0 | 64.3 | 67.5 | **71.2** |
  | spin, post (rad/s) | 0.0 | 50.7 | 81.8 | **86.9** | 65.8 | 60.3 | 54.6 |
  | descent, pre (m/s) | 1.66 | 2.00 | 2.76 | 2.83 | 2.87 | 2.93 | 3.02 |
  | descent, post (m/s) | 1.66 | 13.66 | 13.42 | 13.47 | 9.24 | 6.91 | 6.09 |
  | cone, pre (°) | 51 | 47 | 66 | 68 | 68 | 68 | 68 |
  | cone, post (°) | 51 | 88 | 87 | 87 | 78 | 72 | 73 |

  Two readings. **The twist mechanism does not depend on the cut terms** — a centred-CoM twisted
  seed still spins up on its own, so shape-driven autorotation was never an artifact. But the
  smooth monotonic spin-up (pre-cut: every case `spiral`, all but one converged) becomes a
  non-monotonic curve peaking at 15°, and **descent rises 4–5× as the seed rolls edge-on**
  (cone 68°→87°). These are converged attractors at 5–15°, not integrator noise.

  This is the cleanup's central measurement, and it is a useful one: it quantifies how much of
  the model's plausible-looking behaviour was scaffolding. A phantom moment worth 25% of the
  torque budget, plus roll damping running at 2× reality, were together holding the cone
  shallow. With both gone and nothing physical yet replacing the lift they stood in for, the
  seed has no mechanism to resist going edge-on. **The 3D model is not usable for quantitative
  work until phase 4 lands** — that is expected, and it is the argument for phase 4 rather than
  against the cuts, which removed provably-wrong terms.
- [ ] **3 — Benchmark gate: what survives?** *(run — confirm before ticking.)* Snapshot
  `model_test_results/2026-09-20_134952_af98e91`, now with a fourth stage, `shape3d/`, that
  sweeps twist and curvature (four families × 0–35° tip angle, centred nut) alongside the three
  planar stages. `generateModelBaseline` records both models' switch sets, since they now differ.

  **Planar is frozen — verified three ways.** Against the five-week-old `a51df23` baseline, the
  400-cell mode grid is byte-identical, the test suite summary is byte-identical, and all three
  CoM-movement scenarios reproduce their dwell-mode sequences exactly.

  **The span force is measurably inert.** The flat reference (span force ON) and the 0° curvature
  rows (span force OFF) are otherwise the same seed. They differ by ≤ `7e-4` relative on every
  metric — descent `1.6589` vs `1.6588` m/s, cone `50.95°` vs `50.99°`. It can go in phase 4
  without ceremony.

  **The failure is an ATTITUDE failure, not a lift-magnitude one.** The collapsed cases reach
  genuine terminal equilibrium (vertical aero force / weight = `1.00`) while descending 12–14 m/s
  at cone 84–89°, i.e. edge-on — the plate normal is roughly horizontal, so the large sectional
  forces point sideways and cancel around each revolution, leaving almost nothing supporting the
  weight. Time-weighted on a uniform grid, the strip angle-of-attack distribution is **nearly
  identical across every case** — attached ~10%, LEV band (20–60°) ~22%, bluff (>75°) ~57% —
  including the flat reference descending at 1.66 m/s. The sections are not operating
  differently; the seed is *oriented* differently.

  **And the cuts mostly made the model BETTER.** Running all three shared stages under both
  models on the same flat seed (`mode_grid/` vs `mode_grid_3d/`, etc.) and diffing with
  `compareBaselineModels`:

  | | planar | shape3d | delta |
  |---|---|---|---|
  | diving cells | 88 | 56 | **−32** |
  | autorotation cells | 163 | 178 | +15 |
  | spiral / tightSpiral | 83 / 55 | 95 / 62 | +12 / +7 |
  | mean descent (unchanged cells) | 5.27 m/s | 4.20 m/s | **−1.07** |
  | mean cone | 66.3° | 61.2° | **−5.1°** |

  27.5% of the 400 cells changed mode, and the dominant migration is `diving → autorotation`
  (×33) and `diving → spiral` (×22). Slower descent, shallower cone, less diving, more
  autorotation — every shift is *toward* real samara behaviour. The test-suite summaries are
  identical (though that file records only status/descended, so it is a coarse check), and of the
  three CoM-movement scenarios only the spanwise sweep changed, again gaining rotation
  (`parachuting → diving → autorotation → tightSpiral → autorotation`).

  **The regression that matters: which regime breaks.** Of the seven documented mode-eliciting
  inputs, planar scores 7/7 and shape3d **6/7**. The single failure is `flutter+spiral`, whose nut
  sits at `0.01·S` — a 0.5 mm spanwise offset, i.e. an essentially **centred CoM** — and it
  collapses edge-on to 13.65 m/s at 89.2°, exactly like the centred-CoM twist cases. Every case
  with a substantial CoM offset is fine or improved.

  That is the diagnosis: **a seed with an offset CoM gets its restoring moment from the real
  weight × arm couple, and the cuts helped it. A near-centred or shape-driven seed had no such
  couple, and the phantom span torque was standing in for the missing one.** Which is precisely
  the moment that spanwise load redistribution supplies — see phase 4.
- [ ] **4 — Then add real physics back.** Phase 3 **reordered this list.** The failure mode is
  attitude equilibrium, not sectional lift magnitude, so the items that act on the *moment*
  balance now come first:
  1. **Induced inflow + tip loss.** The strongest candidate, and phase 3 pinned down exactly why.
     What the cuts removed was a pitching/rolling moment (the span torque was 25% of the moment
     budget; `Tx` halved the roll damping), and the cone equilibrium went with it — but *only for
     seeds with a near-centred CoM*, which have no weight × arm couple to supply that moment
     instead. Spanwise load redistribution is the real mechanism that does, and it is the honest
     version of the span-CoP migration that was cut: it shifts the CoP inboard steadily rather
     than oscillating at spin frequency. AR = 3.33, so the elliptic slope is 67% of the 2D value.
     **Success criterion:** `flutter+spiral` (nut at `0.01·S`) recovers `spiral` instead of
     collapsing edge-on, and the centred-CoM twist cases stop diving at 13 m/s — without
     undoing the mode-grid gains on offset-CoM cells.
  2. **LEV augmentation.** Still real and still needed (the separated branch caps at
     `CT = 0.95`; LEV-bearing samara sections reach ~2), but phase 3 argues it will not fix the
     descent on its own: multiplying a force that is pointing sideways and cancelling does not
     restore vertical support. It changes the torques too, so the two are coupled — implement
     inflow/tip-loss first, then judge LEV against a model that holds a sensible cone.
  3. **`-Ȧv`.** Free, first-principles, 0.77 g. Deferred out of phase 2 to keep the benchmark
     attributable; no reason left to defer it further.
  4. **Nut form drag.** Real but least relevant to this failure — for a centred nut it acts near
     the CoM and contributes little moment.
  5. **Edge crossflow drag**, replacing `computeSpanForce`, which phase 3 measured as inert.
  LEV and the aspect-ratio correction push in opposite directions; both are real and they do not
  cancel.

### Bugs — not working to spec, or broken in a meaningful way

- [ ] `Tx` roll damping double-counts the strips' own roll moment — **now proven, not suspected.**
  A pure-roll state (ω_x only, v = 0, CoM centred) isolates the strips' own roll moment; it
  converges on `Tx` as the discretisation refines: ratio `1.0141` at 12 strips → `1.0009` at 48
  → `1.0001` at 192, and holds at `1.0002` across ω_x = 5→60. Both reduce to the same closed
  form, `-ρ·c·CD₂·ω|ω|·S⁴/64`. Roll damping is currently **2× reality**. `spanSpinDamping`'s own
  docstring spots the issue (span *is* the axis the strip loop already discretises) and draws the
  wrong conclusion — call it once, rather than not at all. Resolution: cut it (phase 2 above).
- [ ] The span **torque** is the strips' own normal force re-applied at a fictitious span CoP.
  `seed6DOFODE` discards `F_span_full(2)` from the force sum as a double count, then crosses it
  with a ~15 mm arm anyway. Measured over a settled autorotation: `tau_span` is **25.1%** of
  `|tau_body|`, and it splits `6.861e-05` from the discarded y-force against `1.070e-07` from the
  force actually applied — i.e. **100% of it comes from a force that is never applied**, by a
  factor of 640. Meanwhile the force it exists to supply is `1.641e-05` N, **1.0% of weight**.
  A moment with no corresponding force matches no pressure distribution. Resolution: cut (phase 2).
  This also falsifies the justification comment in `seed6DOFODE.m` — strip theory *can* produce
  roll from span-offset loading, since `cross(r_cp, F)` has an x-component whenever `zGeo ≠ c_z`.
- [ ] `translationDynamics` drops the added-mass rate term `-Ȧv`, whose docstring calls it "tiny
  for a seed in air". It is not: finite-difference-verified to `6.2e-09` against the analytic
  form `ω×(Av) - A(ω×v)`, it runs `1.389e-03` N in settled autorotation — **31% of the net force,
  82% of weight, 0.77 g of acceleration error**. Added mass is only 6.3% of `M`, but the *rate*
  term is amplified by ω (~19 rad/s). Cheap to add and first-principles; put it behind a switch
  so the planar baseline stays reproducible.
- [ ] `normalSpinDamping` (`Ty`) is dimensionally a **force**, not a torque: `ρ·R⁴·ω²` is N, and
  `Tr`/`Tx` each carry a multiplicative width (`dz`, `chord`) that `Ty` lacks. That makes it ~20×
  an honest skin-friction estimate — but it contributes 0.1% of the torque budget, so it is
  harmless and does nothing. Resolution: cut with the rest (phase 2) rather than fix.
- [ ] `testing/planar/runSeedModeGrid.m` hardcodes `seed6DOFODE` instead of dispatching through
  `seedRHS(sp)`, so it silently ignores `cfg.shapeModel` and cannot run the 3D model at all.
  Blocks the phase-3 benchmark; fix first. (It is also serial — worth a `parfor` at 40×40.)
- [x] `classifyFlightMode` thresholds calibrated — all seven modes now label correctly
  against the known working inputs (autorotation vs. spiral split on cone steadiness).
  Caveat: verified only against the txt-file inputs; not tested across a broad range of
  autorotations, and may need re-calibration for future physics (e.g. LEV lift).
- [ ] Autorotation descends too fast and at too steep a cone vs. real samaras (blade-element
  theory under-predicts lift; see LEV feature below).
- [x] `testing/Working Dynamics Inputs 7-24-26.txt` inline comments said "(false by default)"
  for enablers that actually default `true` — fixed to match the code (`enableSpanForce`,
  `enableSpanCOPMigration`, `enableTxDamping` now read "(true by default)"; the two that really
  default false were left as-is). The newer 8-2-26 inputs file was already accurate.

### Features — still to implement

- [ ] Leading-edge-vortex (LEV) lift augmentation — a sectional Polhamus-style term to fix
  the autorotation lift deficit (the physics strip theory omits). The separated branch caps at
  `CT = 0.95`; real LEV-bearing samara sections reach ~2. *Phase 4, highest value.*
- [ ] Finite-aspect-ratio / tip-loss correction. AR = S/c̄ = 3.33, so the elliptic lift slope is
  `a₀/(1+a₀/πAR)` = **67% of the 2D value** — currently uncorrected. This is also the honest
  version of the span-CoP migration being cut: real spanwise load redistribution shifts the CoP
  inboard *steadily*, it does not oscillate at the spin frequency. *Phase 4.*
- [ ] Nut form drag — the nut is currently a drag-free point mass, though in autorotation it sits
  off-axis and sweeps at ~19 rad/s. A real force and moment that are entirely absent. *Phase 4.*
- [ ] Edge crossflow drag as the honest replacement for `computeSpanForce`, *if* phase 3 shows the
  seed needs any spanwise resistance: `F_z = -½ρ·C_d·(t·c̄)·|v_z|v_z` at the area centroid — one
  constant, pure drag, no borrowed lift or CoP, moment from the genuine CoM arm. The current law
  has the wrong shape as well as the wrong size: at β = 0 (pure spanwise slide, the gap it was
  written to fill) it uses the attached `CD0` and returns *less* than a naive edge estimate, while
  at β = 45° it returns 4× that estimate — largest where the strips already cover it, smallest
  where it is needed.
- [ ] 2D validation suite vs. `minimal_imp`: RHS-match at matched states (integration-free)
  plus trajectory comparison, to lock the baseline before further tuning.
- [x] Calibrate `classifyFlightMode` thresholds against the known mode-eliciting inputs.
  (Done for the txt-file inputs only; not validated across a broad range or for future
  physics such as LEV lift.) Some follow-up tuning has been done to smooth the phase map:
  a `tiltChaos` gate so a run that merely fails to settle is no longer forced to `chaotic`
  (only genuinely incoherent tumbling is), a longer default settling time (`tspan` 12 s),
  and a helix-radius fix so near-straight tracks no longer report a spurious finite radius.
  This cut the phase-map "islands" from ~51 to ~31 (the rest are real spin-onset bistability).
- [x] Flight-mode phase map: sweep a chord×span grid of nut positions, classify each, and
  render a 2D mode diagram (`runSeedModeGrid`). Caveat: it's an in-plane (chord, span) slice
  at modest resolution — parachute (out-of-plane nut) doesn't appear, and the large
  autorotation region reflects the current lift deficit, not the true samara map.
- [ ] Standalone intermediates-recovery helper (rebuild per-strip / whole-seed quantities
  over a saved trajectory, since `ode45` discards the second output).

### Outputs — what the code needs to produce, and whether it works yet

- [x] 3D CoM trajectory + Euler-angle history (`visualizeSeedTrajectory`).
- [x] Frame-by-frame seed animation with local velocities (`animateSeed`).
- [x] Per-case and per-group sweep figures (`runSeedTestSuite`).
- [x] Trajectory-metrics extraction (`computeTrajectoryMetrics`).
- [x] Torque-budget + span-CoP/CoM/`F_span` diagnostic animation (`runSpanForceComparison`).
- [x] Time-varying-CoM (moving-nut) test with commanded-path plots (`runComMovementTest`).
- [x] Automatic flight-mode classification — calibrated and since tuned (added a `tiltChaos`
  gate so non-convergence alone no longer reads as `chaotic`, plus a longer settling time and
  a helix-radius fix) to smooth the phase map; labels all seven known modes correctly.
  *Caveat: verified against the txt-file inputs only, not a broad range of autorotations, and
  may need re-calibration for future physics (e.g. LEV).*
- [x] Flight-mode phase diagram (chord×span map) via `runSeedModeGrid`. *Caveat: in-plane
  slice, modest resolution; autorotation region inflated by the lift deficit.*
- [ ] Compiled PDF of `derivations/seed6DOF_physics.tex` for inline viewing on GitHub.

### Extra — optional, future

- [ ] Tune the span torque (`enableSpanTorque` / `aero.C_span_torque`) — **likely superseded.** The
  audit showed the span torque is a phantom moment (see *Bugs*), so it is being cut in phase 2
  rather than tuned. The trade it was balancing (fluttering+spiral vs. autorotation steepness)
  goes to phase 3/4: measure what survives the cut, then resolve it with LEV and tip loss.
- [ ] Prune the now-dead aero constants (`C_span_torque`, `k0_spanTorque`, `C_Tx`, `C_fy`) and the
  `spanSpinDamping` / `normalSpinDamping` files once the planar model is retired. They stay on
  disk through phases 2–3 so the planar baseline keeps running bit-identically.
- [ ] *Not worth doing:* the `(π-|α|)` branch quirks preserved verbatim from the 2D code cost at
  most **0.7% of peak CT and 0.2% of peak CD** over the full ±180° range (the tanh blend
  suppresses them). Measured, not estimated — keep the code comments, skip the fix.
- [ ] Momentum-rigorous internal mass movement (the moving-nut test is a kinematic
  prescription only — no reaction of the sliding mass).
- [ ] Tapered / arbitrary planform support beyond the rectangular test seed.
- [ ] Non-planar (3D) seed shape — **built, not yet experimentally validated.** Lives in its
  own `physics3d/` module (`seedParams.model = 'shape3d'`) beside the frozen planar model;
  both feed the shared rigid-body core `physics/rigidBody6DOF.m` and are selected by
  `cfg.shapeModel` + `testing/helpers/seedRHS.m`. What works:
  - **Builder** `physics3d/setupSeedShape3D.m` — spanwise **twist** θ(s) (pitch about the span
    axis) and **curvature** φ(s) (dihedral about the chord axis), which compose. The flat span
    coordinate is treated as ARC LENGTH and the tangent integrated (`dz=cosφ ds`, `dy=sinφ ds`)
    to place each strip's physical `(z,y)`, so strip widths stay arc-length-correct.
  - **Mass/inertia** — curvature-correct: the wing centroid gains a `y` shift and each strip's
    local inertia is rotated into its own frame before the parallel-axis sum. Verified
    analytic-exact against the closed-form thin plate (~1e-14) and shown to reduce to the
    planar builder.
  - **RHS** `physics3d/seed6DOFODE3D.m` — projects each strip's velocity into its local
    (chord, normal) frame and rotates the force back to body; `computeSeedLocalVel` now
    transports to the strip's true 3D position.
  - **Toggles** — the planar `computeSpanForce` hack auto-disables for a *curved* seed (it
    would double-count the spanwise force the tilted strips now produce from geometry);
    `enableAddedMass3D` switches on the full per-strip added-mass tensor
    `A_trans = Σ mᵢ(n̂ᵢn̂ᵢᵀ)` (off by default; reduces exactly to the flat form).
  - **Viz + tests** — `visualizeSeedShape` draws per-strip quads (shows twist, curvature and a
    varying chord); `testing/shape3d/` holds the twist and curvature suites plus two
    regressions: flat equivalence (bit-identical to the planar model) and the curved
    mass/inertia checks.
  - Twist reproduces the observed effect: a **centred-CoM** twisted seed spins up on its own
    (0→71 rad/s over 0–30° tip twist), i.e. shape-driven autorotation with no mass offset.
  - Curvature behaves as the *complementary* mechanism, and the suite confirms the split: a
    **symmetric bowl produces no spin at all** (0.00 rad/s across 0–35° tip dihedral) but
    collapses the cone 51°→0.1° and turns the seed into a parachute, while **asymmetric**
    (single-side) curvature does spin it. The out-of-plane CoM shift is curvature's signature —
    +1.40 mm at 35°, matching the independent hand calculation in `testCurvedMassProperties`.
  Remaining: validate against real curved/twisted seeds (none of this 3D aero is experimentally
  checked), and optionally accept a general YZ polyline profile instead of the current
  parametric twist/curvature functions.
- [ ] Unsteady CoP lag (first-order relaxation as a real state) instead of the stateless
  reduced-frequency attenuation.
- [ ] Crossflow-drag-only span-force variant (pure drag, no borrowed lift / CoP migration).
- [ ] Full 2D spanwise re-discretization (the "combined" two-plane model, previously set aside).
