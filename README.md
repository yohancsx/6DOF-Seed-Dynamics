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
├── Seed_Dynamics_ODE_Test.mlx     ← START HERE: live script that runs a simple drop
├── physics/                       core dynamics
│   ├── seed6DOFODE.m              main ODE right-hand side (the integrator function)
│   ├── setupSeedShapeAndMass.m    builds the seed: strip geometry, CoM(t), inertia(t)
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
├── visualization/                 visualizeSeedTrajectory.m, animateSeed.m, ...
├── testing/                       test-suite scripts (see below) + testing/helpers/
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

> **Note on the spanwise-flow additions.** `computeSpanForce`, `spanSpinDamping` (`Tx`),
> and `normalSpinDamping` (`Ty`) are experimental whole-seed extensions gated behind
> `seedParams.enable*` switches. The validated, 2D-reducible baseline runs with the
> spanwise force off; see the derivation for what each switch does.

---

## Running the code

**Requirements:** MATLAB (R2020b or newer recommended). No toolboxes required beyond core
MATLAB. But it's nice if you have the aerospace toolbox for quaternion conversion.

### Quick start — a single drop

Open and run **[`Seed_Dynamics_ODE_Test.mlx`](Seed_Dynamics_ODE_Test.mlx)**. It:

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

| Script | What it does |
|---|---|
| `runSeedTestSuite.m` | Parameter sweeps (nut mass/position, initial tilt/spin, asymmetry, strip-count convergence); saves a trajectory figure per case and a per-group overlay. |
| `runSeedModeSuite.m` | Elicits and auto-classifies the biological flight modes (glide, dive, spiral, autorotation, parachute). |
| `runSpanForceComparison.m` | Four hand-tuned cases isolating the spanwise-force physics, with a torque-budget diagnostic. |
| `runComMovementTest.m` | Drives a **time-varying** CoM (the nut slides within the body) to test mode transitions. |

Each script has an editable configuration block at the top and assumes `physics/` and
`visualization/` are on the MATLAB path.

---

## Reproducing the flight modes

The configuration below reproduces the full set of biological descent modes. **It is now
the code default** ("FULL-minus-geomVelocity"), set in `setupSeedShapeAndMass` /
`buildSeedParams` / `seed6DOFODE`, so a seed built the normal way already uses it — you do
not need to set these by hand. Recorded in
[`testing/Working Dynamics Inputs 8-2-26.txt`](testing/Working%20Dynamics%20Inputs%208-2-26.txt).

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

### Bugs — not working to spec, or broken in a meaningful way

- [ ] `Tx` roll damping double-counts the strips' own roll moment (strips offset in span
  already produce roll damping via their normal-load distribution). It's on in the working
  config; decide partial-strength vs. off vs. strip-only, and settle the derivation.
- [x] `classifyFlightMode` thresholds calibrated — all seven modes now label correctly
  against the known working inputs (autorotation vs. spiral split on cone steadiness).
  Caveat: verified only against the txt-file inputs; not tested across a broad range of
  autorotations, and may need re-calibration for future physics (e.g. LEV lift).
- [ ] Autorotation descends too fast and at too steep a cone vs. real samaras (blade-element
  theory under-predicts lift; see LEV feature below).
- [ ] `testing/Working Dynamics Inputs 7-24-26.txt` inline comments say "(false by default)"
  for enablers that actually default `true` in the current code — misleading, should match.

### Features — still to implement

- [ ] Leading-edge-vortex (LEV) lift augmentation — a sectional Polhamus-style term to fix
  the autorotation lift deficit (the physics strip theory omits).
- [ ] 2D validation suite vs. `minimal_imp`: RHS-match at matched states (integration-free)
  plus trajectory comparison, to lock the baseline before further tuning.
- [x] Calibrate `classifyFlightMode` thresholds against the known mode-eliciting inputs.
  (Done for the txt-file inputs only; not validated across a broad range or for future
  physics such as LEV lift.)
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
- [x] Automatic flight-mode classification — calibrated; labels all seven known modes
  correctly. *Caveat: verified against the txt-file inputs only, not a broad range of
  autorotations, and may need re-calibration for future physics (e.g. LEV).*
- [x] Flight-mode phase diagram (chord×span map) via `runSeedModeGrid`. *Caveat: in-plane
  slice, modest resolution; autorotation region inflated by the lift deficit.*
- [ ] Compiled PDF of `derivations/seed6DOF_physics.tex` for inline viewing on GitHub.

### Extra — optional, future

- [ ] Tune the span torque (`enableSpanTorque` / `aero.C_span_torque`): it enables the
  gentle fluttering+spiral mode but steepens and speeds up autorotation. Find a middle
  setting, or resolve the trade with LEV lift so both modes coexist.
- [ ] Momentum-rigorous internal mass movement (the moving-nut test is a kinematic
  prescription only — no reaction of the sliding mass).
- [ ] Tapered / arbitrary planform support beyond the rectangular test seed.
- [ ] Non-planar (3D) seed shape: define a profile as segments in the body **YZ** plane
  (normal × span) and extrude it by a chord width along **X**, giving the wing real
  dihedral/camber/curl. Relevant for actual samaras (they are not flat), but can wait until
  the larger flat-plate results are in. The rigid-body core is already 3D; the flat-plate
  assumption is confined to strip-level geometry and aero. Each strip's local frame is the
  body frame *rolled about the chord (X) axis* by its local dihedral angle φ, so the change
  is: give each strip a `y` position + roll angle, rotate velocities in / forces out about X.
  General places to touch:
  - **New builder** (the shape-gen script): from the YZ polyline + chord, emit per-strip
    `xgc, ygc, zgc`, arc-length width `ds` (not the z-projection), and roll angle φ / normal
    `n̂`, into the same `strips.*` struct the ODE already consumes (flat = the φ=0 case).
  - `physics/setupSeedShapeAndMass.m` — plate centroids get a real `y` (lines ~184, ~244);
    roll each plate's local inertia tensor by φ (`I_local → Rᵢ·diag(...)·Rᵢᵀ`). Parallel-axis
    sum is already 3D.
  - `physics/seed6DOFODE.m` — put `ygc` into `r_cp`/`r_geoCenter` (lines ~303–304); project
    strip velocity onto the local (chord, normal) via φ instead of `vStrip(1),(2)` (~292–293);
    rotate the returned strip force back to body about X; project spin onto the *local* span
    axis for the rotational-lift `omega_z`.
  - `physics/helpers/computeSeedLocalVel.m` — `0 → ygc(i)` (line ~97).
  - `physics/mass/getAddedMass.m` — broadside added mass acts along each `n̂ᵢ`, so
    `A_trans = Σ mᵢ·(n̂ᵢ n̂ᵢᵀ)` (full 3×3) and `A_rot` needs the y-lever arms; the flat
    approximation is a small error for gentle curl and can be deferred.
  - `visualization/visualizeSeedShape.m` / `animateSeed.m` — render the extruded curved sheet
    instead of the flat x-z plate.
  - Design note: a genuine 3D shape produces the spanwise force from geometry, so the ad-hoc
    `computeSpanForce` / span-torque block becomes redundant — default `enableSpanForce = false`
    for 3D seeds. `computeStripForces` / `computeAngleOfAttack` need no change if the ODE does
    the projection (treat their output as the strip-local frame).
- [ ] Unsteady CoP lag (first-order relaxation as a real state) instead of the stateless
  reduced-frequency attenuation.
- [ ] Crossflow-drag-only span-force variant (pure drag, no borrowed lift / CoP migration).
- [ ] Full 2D spanwise re-discretization (the "combined" two-plane model, previously set aside).
