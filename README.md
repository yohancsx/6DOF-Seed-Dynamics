# 6-DOF Seed Dynamics

A MATLAB simulation of the free flight of an autorotating winged seed (a **samara**) —
modelled as a flat plate with a discrete "nut" mass — using **quasi-steady strip-theory
aerodynamics**. It is a full 6-degree-of-freedom (13-state) extension of the classic 2D
Andersen–Pesavento–Wang falling-plate model into three dimensions.

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
| **[`Olivia Code/minimal_imp_Commented.m`](Olivia%20Code/minimal_imp_Commented.m)** | The **2D reference model** (Andersen–Pesavento–Wang falling plate). Ground truth for all aerodynamic physics — every coefficient and force term in the 3D code traces back to it. |
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
MATLAB.

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

## The model in brief

- **Per strip:** translational lift + drag (at the migrating centre of pressure) and
  rotational lift (at mid-chord), from the 2D quasi-steady coefficient laws.
- **Whole seed:** gravity, added mass, gyroscopic and time-varying-inertia torques, and the
  optional spanwise-flow force/torque and spin-damping terms.
- **Deliberately omitted:** buoyancy (negligible in air) and added-mass CoM-offset coupling.

Full derivation, including every coefficient branch and added-mass term, is in
[`derivations/seed6DOF_physics.tex`](derivations/seed6DOF_physics.tex).

---

## References

The model and its validation targets draw on:

1. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Unsteady aerodynamics of fluttering
   and tumbling plates.* Journal of Fluid Mechanics, 541, 65–90. — the 2D quasi-steady
   falling-plate model this code extends.
2. Andersen, A., Pesavento, U., & Wang, Z. J. (2005). *Analysis of transitions between
   fluttering, tumbling and steady descent of falling cards.* Journal of Fluid Mechanics,
   541, 91–104.
3. Pomerenk, O., & Ristroph, L. (2024). *Aerodynamic equilibria and flight stability of
   plates at intermediate Reynolds numbers.* Journal of Fluid Mechanics.
   [arXiv:2408.08864](https://arxiv.org/abs/2408.08864) — steady equilibria (gliding,
   diving) and stability of the falling-plate modes.
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
