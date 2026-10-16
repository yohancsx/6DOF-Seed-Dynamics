# Handoff: import the 29 measured seeds and add a real-seed test stage

This is a task brief for a separate Claude instance. It holds the data facts, the plan,
and every place the work connects to the current code. Work through it in the order
given in §6.

**Nothing here is decided until the user says so.** §5 lists the open decisions, each
with a recommended default. Confirm them with the user before writing code.

> **UPDATE 2026-10-15: a cleanup landed after this brief was written. Read this first.**
> The sections below have been updated to match, but in short:
> - **One model.** The planar model is archived (`archive/planar/`, git tag
>   `planar-final`). `buildSeedParams` is shape3d-only (`cfg.shapeModel` defaults to
>   `'shape3d'`; `'planar'` is an error), and so are `seedRHS` and `validateSeedParams`.
> - **The builder no longer wraps the planar builder.** `setupSeedShape3D` has four
>   sections:
>   1. STRIPS: `physics3d/slicePlanform.m` gives the flat strips and flat mass;
>   2. FRAMES: twist and curvature;
>   3. MASS: `curvedMassProperties` when curved;
>   4. PHYSICS: switches, LEV, `seedThickness`, `.model`, `.baseSeedParams`.
>
>   **The `stripGeom` path replaces sections 1–3. Section 4 is the shared tail, already
>   factored**; there is no dead-switch stripping any more.
> - **`runSingleMode` now labels with the six-part classifier** (`result.mode`,
>   `result.parts`, `result.metrics` from `computeSeedModeMetrics`). The old
>   `classifyFlightMode` is archived; don't use it.
> - **Tracking:** use `testing/helpers/appendTrackingRow.m`. It merges new columns into
>   `tracking_history_3D.csv`, so decision 6 below is settled: add `real_*` columns to the
>   main history.
> - **Tests:** every regression test is a function returning `res` (`.name/.pass/.detail`).
>   `testing/runRegressionTests.m` runs them all, and it is stage 0 of the suite (warn
>   only). Add `testRealSeedLoader` to its list. `testMinimalImpEquivalence` is the primary
>   flat test.
> - **Guard rails:**
>   - `testing/shape3d/testRefactorInvariance.m` pins the model bit for bit. Run
>     `('compare')` after the builder change: the existing seeds must stay 44/44 identical.
>   - `testActivePathClean` fails if active code calls archived code.
> - **Paths:** scripts find the repo through `seedPaths.m` (the 3-line bootstrap at the top
>   of every script). `P.seedData` points at the real-seed `.mat`.
> - **Suite:** stages are now 0 and A–F, with `cfg.profile = 'quick' | 'full'`. Stage G
>   goes between F and 8b as planned. The grid picker moved to `testing/tools/`.

---

## 0. Ground rules for this repo

These are standing user rules. Follow them.

- **Free-flight data is validation only.** Never tune a coefficient to make a real seed
  match its trials. That covers the user's seed trials and Hou et al. 2025.
  - A new physics term needs a form from theory and a magnitude from independent static
    data.
  - This stage compares the model with real seeds. It does not calibrate anything.
- **Use the least "fake" physics.** Don't add new aero terms to make real seeds behave.
  Missing physics, such as the kernel's own aerodynamic body, gets documented as a known
  limitation.
- **Keep generated binaries out of git** (videos, figures, `.mat` results).
  - `model_test_results/` is gitignored except for its README.
  - The 24 MB seed `.mat` stays outside the repo. Reference it by path.
- **README roadmap:** after finishing an item, ask the user before marking it `[x]`.
- **The user commits.** Don't commit unless asked. Keep the tree clean before any
  baseline run, because run folders are stamped with the git hash plus `-dirty`.
- **Running MATLAB:** `matlab -batch "run('<script>.m')"`. Inline quoting through
  PowerShell breaks easily, so put multi-line work in a scratch `.m` file in the session
  scratchpad and run that file.
- **Code style:**
  - Match the existing header-comment style: a long `%` header with INPUT, OUTPUT and the
    reasoning.
  - Body axes are **x = chord, y = normal, z = span**.
  - World **Y is vertical**.
  - Quaternions are scalar-first, body→world.

---

## 1. The data

- **File:**
  `C:\Users\yohan\OneDrive\Documents\Research Stuff\Seed Dynamics Code\Real Seed Validation Data\seedResults.mat`
- **Variables:**
  - `R`: a 1×29 struct array, one entry per seed.
  - `T`: a 29×54 summary table.
  - `opts`: the fitting options. It sets wing density `0.2581` g/cm³ and kernel density
    `1.107` g/cm³. These are assumed values, not measured per seed.
- **Units are mm, g and g·mm².** Convert to SI:
  - lengths ×1e-3;
  - mass ×1e-3;
  - inertia ×1e-9.
- **Seeds 3–31 are present** (no seed 1 or 2). Names look like `seed10MeshCleaned`.
  `R(k).seedNum` holds the number.

### Per seed `R(k)`

| field | contents |
|---|---|
| `.align` | `R` (3×3 alignment rotation), `c0`, `extents`, `bboxMin/Max`. Check `det(R) = +1`; a reflection would mirror the seed's chirality. |
| `.metrics` | `length`, `width`, `planformArea`, `meanChord`, `aspectAero`, `thickMeanWing`, `curvatureSagitta`, `twistRange`, `twistRatePerSpan`, `twistLinearR2`, `asym.*`, `half.*`, … |
| `.plates` | **the fitted shape**, described below |
| `.lengthCheck` | `measuredMm`, `meshMm`, `errPct` |
| `.figures` | paths to the user's fitting figures (alignment, sections, profiles, plates). Useful for eyeballing the import. |

### `R(k).plates`

- **`numPlates = 12`**, `plateShape = 'trapezoid'`.
- **`hinges`** (13×3): the hinge points between plates.
- **`plates(i)`, a 12×1 struct.** Each plate is a FLAT trapezoid with these fields:
  - **Position and size:** `x0`, `x1` (span limits); `H0`, `H1` (hinge points); `length`
    (along `e1`); `chordIn`, `chordOut`; `centroid` (1×3); `corners` (4×3).
  - **Frame:** `e1` (span axis), `e2` (chord axis), `n` (normal), all unit vectors.
    I checked plate 6: `e2 × n = e1`, so the frame is right-handed.
  - **Shape angles (deg):** `twistDeg`, `dihedralDeg`, `sweepDeg`.
  - **Mass properties:** `thickness`, `area`, `areaProj`, `volume`, `mass` (g).
  - **`inertia`:** a 3×3 tensor about the plate's own centroid, in **mesh axes**, g·mm².
    I verified it equals density × (sum of the two other second moments). Use it as an
    **independent cross-check**.
  - **`Jlocal`, `Jaligned`:** VOLUME second moments (∫x² dV), not inertia. Don't use them
    as inertia.
- **`.joints`:** `position`, `axis`, `foldDeg`, `twistJumpDeg`, `splayGap`.
- **`.kernel`:**
  - `volume`, `centroid` (1×3, mm), `centroidS`, `sRange` (in normalised span `s ∈
    [-1,1]`);
  - `sectionArea` (61×1, the kernel cross-section area at each station);
  - `region` (61×61 mask), `equivalentRadius`, `density`, `mass` (g).
- **`.fit`:**
  - `modelMass` = sum of the plate masses + the kernel mass;
  - `measuredMass`, `massErrPct`;
  - `volumeRatio` (≈1.17: the middle plates double-count some kernel volume);
  - residuals.
- **`.ref`:** 61-station reference profiles (`x`, `yRef`, `zRef`, `twistMid`, `chord`,
  `wingThickness`). It also holds function handles (`fy`, `fz`, …) that depend on a
  `holdInterp` function that may not be on the path, so **don't rely on the handles**.
- **Station grid:** `R(k).sections.xStation` (61×1, mm) maps `s` to span position, with
  `halfSpan` alongside.

### Mesh axes

- **Mesh frame:** x = span (length), y = chord (width), z = thickness/normal.
- **To body:** `body = P * mesh` with `P = [0 1 0; 0 0 1; 1 0 0]`, which gives
  `[y; z; x]`.
  - This is a cyclic permutation, so it is a proper rotation.
  - Apply `P` to positions and to the frame vectors, so `chordDir = P*e2`,
    `normalDir = P*n`, `spanDir = P*e1`.
  - Rotate inertia tensors as `P*I*P'`.
- **The sign of +x (span) is arbitrary between seeds** (`metrics.signed.note`, flag
  `T.pcaFlipped`).
  - Dynamics don't care, as long as the seed isn't mirrored.
  - Every seed has negative sagitta, so the bow direction is consistent: tips toward
    +z_mesh, which is body +y.

### What the population looks like (for sanity checks)

- **Size:** span 40–55 mm, chord 14–21 mm, AR 2.6–3.6.
- **Mass:** measured 0.26–0.49 g.
  - The kernel is 0.16–0.36 g, about 75–80% of the total.
  - It sits at mid-span (`|kernelS| ≤ 0.08`) and is spindle-shaped: about 14 mm long
    along the span, with an equivalent radius of about 3 mm.
- **Wing loading is about 0.45–0.55 kg/m².** The suite's synthetic seed is about 0.23,
  so expect roughly 1.4× its descent speed.
- **Bow:** sagitta −0.9 to −6.3 mm, tip dihedral up to ±32°.
- **Twist:**
  - 3–27° range along each seed, mostly **not** linear antisymmetric (median R² ≈ 0.25).
  - Often both tips are washed out the same way. Seed 10, for example, has −7° at one tip
    and −9° at the other.
- **Sweep:** some tip plates reach ±30° (seeds 3, 19, 24, 31).
- **Mass mismatch:** the fitted mass is −27% to +14% off the measured mass. Seeds **10,
  19, 27 and 28** are worse than −24%. Flag them in the output.
- **Leading edge:** the model's centre of pressure (`computeAeroCoeffs` →
  `l_cp_frac`, `computeStripCoP`) flips with the flow direction, so we don't need to
  assign one.

---

## 2. Hooks into the current code

The 3D model's seed contract is documented in `physics/validateSeedParams.m`.

### Seed construction chain

```
bsp (base seed params)
  └─ testing/helpers/buildSeedParams.m (bsp, cfg)        [shape3d only]
       ├─ physics3d/setupSeedShape3D.m(struct('baseSeedParams', bsp))
       │     ├─ 1. STRIPS  physics3d/slicePlanform.m (polyshape → strips + flat mass)
       │     ├─ 2. FRAMES  twist (bsp.twist) / curvature (bsp.curvature) → per-strip frames
       │     ├─ 3. MASS    curvedMassProperties (local fn) when curved
       │     └─ 4. PHYSICS switch defaults + LEV constants (levPlanformConstants) +
       │                   seedThickness, .model = 'shape3d', .baseSeedParams = bsp
       ├─ sets sp.rhoFluid, sp.g
       └─ applies explicit cfg switch overrides (enableEdgeDrag, enableAddedMassRate,
          enableAddedMassMoment, enableAddedMass3D, enableLEV, levApplicationPoint, cfg.lev, cfg.aero)
```

### Fields the RHS consumes

The RHS is `physics3d/seed6DOFODE3D.m`, picked by `testing/helpers/seedRHS.m` from
`sp.model`.

- **`strips`** (1×M unless noted):
  - `chord`, `dz`, `area`, `xgc_body`, `ygc_body`, `zgc_body`, `z_body`, `liftMult`,
    `dragMult`;
  - `chordDir`, `normalDir`, `spanDir` (3×M, body coordinates);
  - `twist`, `dihedral` (for information only).
- **`massParams`:** `tSamples`, `com_t` (3×N), `I_G_t` (3×3×N), `I_G_dot_t`, `M_total`.
  `physics/mass/getMassProperties.m(t, sp)` interpolates these and adds added mass via
  `physics/mass/getAddedMass.m`.
- **Other fields:** `seedThickness` (edge-drag frontal area), `lev`, the switches,
  `rhoFluid`, `g`, and optional `aero`.
- **`getAddedMass`** uses the flat-plate form unless `enableAddedMass3D` is true. The 3D
  form sums per-strip normals and needs `normalDir`, `spanDir` and `ygc_body`.
- **Edge drag** uses a whole-seed span axis `[0;0;1]`, the area centroid, and
  `seedThickness × mean chord`.

### Rendering

The animation and visualisation code reads `sp.baseSeedParams.seedShape` (a polyshape,
drawing x = body z (span), drawing y = body x (chord)):

- `visualization/animateModeTrajectory.m` uses its centroid, vertices and bounding box;
- `visualization/visualizeSeedShape.m` and `animateSeed.m` also use it.

So the importer **must supply a planform polyshape**: the plates projected onto the
(body z, body x) plane, merged with `union`. Then the existing animation code works
unchanged.

### Run, classify and output helpers (all in `testing/helpers/`)

| helper | use |
|---|---|
| `runSingleMode(name, nutPos, q0, omega0, cfg, baseBsp)` | **Overrides** `nutMass_t` with `cfg.nutMass` and `nutPos_t`, so it **can't** be used as-is for real seeds. Integrate directly, as stages E and F do (`seedRHS` + `ode45`, `x0 = [zeros(3,1); q0; zeros(3,1); omega0]`, `odeset('RelTol', cfg.odeRelTol, 'AbsTol', cfg.odeAbsTol)`), or add a small real-seed variant. |
| `computeSeedModeMetrics(t, x, opts)` | Needs `opts.spanLength` (set per seed, in m) and `opts.flipHysteresisTurns`. |
| `classifySeedMode(sm, th)` → `[mode, parts, info]` | The six-part classifier of record, with `th = defaultSeedModeThresholds()`. |
| `seedModePalette()` | Names, codes, colours and the physical flag for the 17 labels. |
| `summarizeModeGrid(labels, absTilt, revolving, refFamily)` | Gives the non-physical fraction and tilt statistics; reusable for a list of seeds (pass `refFamily = []` or strings). |
| `computeTrajectoryMetrics(t, x, metricOpts)` | Gives `descentSpeed`, `verticalSpinMag`, `coneAngleDeg`, … |
| `animateCaseVideo(t, x, sp, name, outDir, aopt)` | Writes a windowed video and never throws. |

`physics/helpers/axisAngleToQuat(axis, angle)` builds the release quaternions.

### The suite: `testing/shape3d/runSeedTestSuite3D.m`

- **§1 config:**
  - `outPath = fullfile(root, "model_test_results")`;
  - stage toggles `cfg.runTwist`, `runCurvature`, `runSweeps`, `runComMove`,
    `runModeGrid`, `runPaperGrid`, plus `cfg.animate`;
  - `cfg.tspan = [0 12]`, `cfg.odeRelTol = 1e-6`, `cfg.odeAbsTol = 1e-8`;
  - `cfg.q0Twist = [1;0;0;0]` (level release); `cfg.omega0`;
  - `cfg.seedModeThresholds`; `cfg.trackingFile`.
  - The user's working copy sets `sw.enableLEV = true` and
    `sw.levApplicationPoint = 'forward'` with label `all_enabled`. **Leave those lines
    alone.**
- **§2:** the run folder `runDir = <outPath>/<timestamp>_<githash>[-dirty]_3D_<label>`.
  It also builds the metric options `so`.
- **§3–§8:** stages A–F. Stage A (lines ~204–234) is the cleanest pattern to copy:
  - a per-case loop, then `classifySeedMode(computeSeedModeMetrics(...))`;
  - optional `animateCaseVideo` into `<runDir>/<stage>_3D/`;
  - a struct saved with `save(...,'cfg','effSwitches')`;
  - lines appended to `summary(end+1)`.
- **§8b tracking:** builds a one-row `table` and writes it with
  `appendTrackingRow(cfg.trackingFile, row)`, which merges new columns (older rows get NaN).
  Add `real_*` columns to that row, with NaN when the stage is off.
- **§9:** manifest and summary.
- The grid picker `testing/tools/pickModeGridRuns.m` is not involved.

### Existing tests to copy and to keep passing

- `testing/runRegressionTests.m` runs all of them (~35 s). **Everything must still pass**
  after any `setupSeedShape3D` change, and `testRefactorInvariance('compare')` must stay
  44/44.
- `testing/shape3d/testShape3DFlatEquivalence.m`: the 3D model on a flat seed must equal
  the (archived) planar model's core.
- `testing/shape3d/testCurvedMassProperties.m`: the pattern for a mass-property
  cross-check.
- `testing/shape3d/testPaperMetrics.m`: 17 synthetic classifier cases.

### Docs to update at the end

Remember the ask-before-ticking rule for the roadmap.

- `README.md`: a test-table row for the new stage and test, and a roadmap line.
- `model_test_results/README.md`: list `real_seeds/`.
- `RESEARCH_LOG.md`: a new phase entry with the first numbers.

---

## 3. Plan

### 3.1 Loader: `testing/helpers/loadRealSeeds.m`

`seeds = loadRealSeeds(file, opts)` returns a struct array with one element per seed.

**Fields:**
- `.name`, `.seedNum`;
- `.bsp`: ready for `buildSeedParams(bsp, cfg)` with `cfg.shapeModel = 'shape3d'`;
- `.meta`, holding:
  - measured and fitted mass, `massErrPct`, the applied mass-scale factor;
  - length, AR, wing loading, twist range, sagitta;
  - `flagMass` (|err| > 20%), `pcaFlipped`, figure paths.

**Options**, with recommended defaults that the user still has to confirm (§5):
- `seeds` (subset, default all);
- `massMode = 'kernel'`;
- `kernelModel = 'distributed'`;
- `subStrips = 1`.

**Steps per seed:**
1. **Convert units** to SI.
2. **Turn plates into strips.** For plate i and each sub-strip j of `subStrips`:
   - Plates are flat, so subdividing along `e1` is exact. Interpolate the chord linearly
     from `chordIn` to `chordOut`, and the centre along the plate's mid-line.
   - Frames: `chordDir = P*e2`, `normalDir = P*n`, `spanDir = P*e1`.
   - Position = sub-strip centroid in body coordinates: `xgc_body`, `ygc_body`,
     `zgc_body` (and `z_body = zgc_body`).
   - `dz = length/subStrips` along `e1`; `area = trapezoid area`;
     `chord = area/dz`, so `sum(chord.*dz)` equals the fitted area exactly.
   - `liftMult = dragMult = ones`.
   - Record `twist = deg2rad(twistDeg)` and `dihedral = deg2rad(dihedralDeg)` for
     information.
3. **Kernel**, `kernelModel = 'distributed'`:
   - Treat it as a stack of discs along the span axis at the stations inside
     `sRange`/`region`.
   - Disc k: area `A_k = sectionArea(k)`, mass ∝ `A_k Δx`, scaled to the kernel mass.
   - Disc centre: `[centroid chord, centroid normal]` taken from `kernel.centroid`, at
     span station `xStation(k)`.
   - Own inertia per disc: `(1/2) m r²` about the span axis and `(1/4) m r²` about the
     other two, with `r = sqrt(A_k/π)`.
   - Sum with the parallel-axis term.
   - `'point'` puts all the kernel mass at `kernel.centroid` with no self-inertia.
4. **Mass mode:**
   - `'kernel'` scales the kernel mass so the total equals `measuredMass`.
   - `'uniform'` scales everything.
   - `'asFitted'` leaves the masses alone.
5. **Mass properties**, constant in time (`tSamples = 0`, `I_G_dot_t = 0`):
   - each plate is a thin trapezoid lamina with mass `plate.mass`;
   - its local inertia is rotated into body axes: `P*plate.inertia*P'`, or an equivalent
     computed from the strip dimensions;
   - add the parallel-axis terms and the kernel bodies;
   - produce `com_t`, `I_G_t`, `M_total`.
6. **Planform polyshape** for the renderers: `bsp.seedShape` is the union of the plate
   corners projected to (drawX = body z, drawY = body x).
7. **Other `bsp` fields:**
   - `seedThickness`: area-weighted plate thickness;
   - `numStrips`: 12 × `subStrips`;
   - `seedDensity`: the wing areal density, needed only by legacy code paths.
8. **Check `det(align.R)`.** If it is −1, warn and skip the seed. Seeds that pass are not
   mirrored.

### 3.2 Builder: explicit-geometry input for `physics3d/setupSeedShape3D.m`

Recommended approach:

- **In sections 1–3:** if `bsp.stripGeom` is present (strips, frames, mass properties and
  polyshape, built by the loader), take the strips and `massParams` from it directly.
  This skips `slicePlanform` and the twist/curvature/`curvedMassProperties` steps.
- **Then section 4 runs unchanged.** It is the shared tail: `seedThickness`, the
  physics-switch defaults, the LEV constants, `.model = 'shape3d'`, and
  `.baseSeedParams = bsp` (which the renderers need).
- **The existing paths must stay bit-identical.** Verify with
  `testShape3DFlatEquivalence` and `testCurvedMassProperties`, and with a short smoke run
  of the suite (stage A on a single twist case: before and after must match).
- **Why this over a separate builder:**
  - `buildSeedParams`, cfg overrides, `seedRHS` and later the grid picker all work
    unchanged;
  - there is no duplicated defaults block to drift.

### 3.3 Runner: `testing/shape3d/runRealSeedSet3D.m`

`out = runRealSeedSet3D(seeds, cfg, outDir)` is a function. It also gets a thin
standalone script cell, or a `if ~nargin` default block, so the user can run it alone.

- **Two releases per seed, both with `omega0 = [0;0;0]`.** Confirm the angle with the
  user (§5).
  - **flat:** `q0 = [1;0;0;0]`, so the plates' best-fit mid-plane is horizontal and body
    y is vertical.
  - **tilted:** `q0 = axisAngleToQuat([0;0;1], deg2rad(cfg.realTiltDeg))`, with a default
    of 10°: a small rotation about the span axis.
- **Integrate** with `buildSeedParams(seed.bsp, cfg)` → `seedRHS(sp)` → `ode45`, using
  `cfg.tspan` and the suite tolerances.
- **Metrics per run:**
  - `classifySeedMode(computeSeedModeMetrics(t, x, so_seed), th)`, with
    `so_seed.spanLength` set to the seed's span (m), plus the parts;
  - `computeTrajectoryMetrics`: descent speed, spin rate, cone angle;
  - span tilt, radius per span, glide ratio (from the `computeSeedModeMetrics`
    fields);
  - time to steady motion.
- **Outputs** go to `<runDir>/real_seeds/`:
  - `real_seeds_summary.csv`: one row per seed × release, with geometry and mass meta
    and the flags;
  - `real_seeds_3D.mat`: results without full trajectories, unless `cfg.saveRealTraj`;
  - `summary_real_seeds.txt`;
  - `real_seeds_overview.png`, with three panels:
    1. a label grid (seed rows × {flat, tilted}) in `seedModePalette` colours;
    2. descent speed vs wing loading;
    3. spin rate vs twist range;
  - `real_seeds_geometry.png`, or one page per seed: a 3D render of the imported plates
    and the kernel. This lets the user check the import against `R(k).figures.plates`.
  - Videos only when `cfg.realSeedVideos` is true (default false), via
    `animateCaseVideo`.
- **Tracked numbers:**
  - non-physical fraction;
  - fraction of seeds whose label differs between flat and tilted (release
    sensitivity);
  - median descent speed and spin rate.

### 3.4 Wire into the suite as stage G

- **New config in §1:**
  - `cfg.runRealSeeds = true`;
  - `cfg.realSeedFile = "<path above>"`;
  - `cfg.realSeedOpts` (massMode, kernelModel, subStrips, seeds);
  - `cfg.realTiltDeg = 10`;
  - `cfg.realSeedVideos = false`.
- **New section between §8 (stage F) and §8b:**
  - call `loadRealSeeds` and then `runRealSeedSet3D(seeds, cfg, fullfile(runDir, 'real_seeds'))`;
  - add summary lines;
  - put the tracked numbers into the §8b row (see the CSV caution in §2).
- **If the `.mat` path is missing**, warn and skip the stage. Don't error: the suite has
  to run on machines without the data.

### 3.5 Tests: `testing/shape3d/testRealSeedLoader.m`

This uses the pass/fail print style of `testPaperMetrics.m`.

1. **Contract:** `validateSeedParams` passes for all 29 seeds after `buildSeedParams`.
2. **Geometry:**
   - `sum(chord.*dz)` matches `plates.fit.areaTotal` (×1e-6) to 1e-6 relative;
   - the span extent matches `metrics.length` to about 1%;
   - every frame is orthonormal and right-handed (`cross(c,n) = s`);
   - `det(align.R) = +1`.
3. **Mass:**
   - `M_total` equals `measuredMass` (kernel mode) or `modelMass` (asFitted);
   - the CoM and `I_G` agree within 1% with an independent sum built from each plate's
     own `inertia` tensor (rotated by `P`) plus the parallel-axis terms plus the kernel.
4. **Flat limit:**
   - a synthetic set of 12 flat rectangular plates with a point kernel at the origin,
     fed through `stripGeom`;
   - it must match today's `setupSeedShape3D` flat rectangle with the nut at the
     origin, same dimensions and `numStrips = 12`;
   - compare strips, `massParams` and one RHS evaluation at a random state, to about
     1e-12.
5. **Resolution:** `subStrips` = 1, 2 and 4 on three seeds (one with high twist, one with
   high sweep, one typical). Labels match, and descent speed moves by less than about 3%.
6. **Mirror:** reflect a seed (z → −z in body, rebuilt with right-handed frames). The
   spin direction flips and the descent speed is the same.

---

## 4. Known limitations to document (don't fix them with new terms)

- **Strip theory doesn't see the kernel's aerodynamic body**, a bump about 3 mm thick.
  Its only aero effect is through edge-drag thickness.
- **The middle plates double-count some kernel volume** (`volumeRatio` ≈ 1.17). The
  kernel mass scaling absorbs this in total mass, but it shifts a little mass from the
  kernel to the middle plates.
- **The flat-plate added-mass form** is an approximation at 32° dihedral. See §5,
  decision 4.
- **Tip plates with large sweep:** strip theory treats each plate as a 2D section normal
  to its own span axis.

---

## 5. Open decisions: ask the user before coding

1. **Mass:** rescale the kernel to the measured mass (recommended) / rescale uniformly /
   leave as fitted.
2. **Kernel:** distributed disc stack (recommended) / point mass.
   - Rough numbers for seed 10: the kernel is about a third of the inertia about the span
     axis (pitch): ≈0.7 against ≈1.4 g·mm² for the wing.
   - A point mass sets that to zero, on the axis the paper comparison flagged as the open
     physics question.
3. **"Small rotation about the span":** a 10° initial angle with zero spin
   (recommended), or a small initial spin rate? What size?
4. **Added mass:** builder default (flat form) for stage G, for consistency with the
   suite, plus a one-off comparison with `enableAddedMass3D = true` (recommended). Or
   switch to 3D for this stage.
5. **Run length:** a fixed 12 s like the suite (recommended for classification), or
   matched to the trial drop height? If matched, what height?
6. **Tracking CSV:** *settled by the cleanup:* add `real_*` columns to
   `tracking_history_3D.csv` via `appendTrackingRow`.
7. **Videos:** off by default (recommended), or on?

---

## 6. Order of work

Stop after each step that changes behaviour, report, and let the user commit.

1. **`loadRealSeeds.m`, the `setupSeedShape3D` `stripGeom` path and
   `testRealSeedLoader.m`.**
   - Run the new test, `runRegressionTests` and `testRefactorInvariance('compare')`, and
     report the results.
2. **Geometry check sheet.** The user compares the 29 imports with their fitting
   figures.
3. **`runRealSeedSet3D.m` standalone,** on 3 seeds first and then all 29.
4. **Stage G in `runSeedTestSuite3D.m`** and the tracking hook.
   - Smoke-run with the other stages off and `cfg.tspan` short, into a scratch
     `outPath`.
5. **Docs** (README test table and roadmap line, `model_test_results/README.md`,
   `RESEARCH_LOG.md`), and ask before ticking anything.
6. **Baseline:** once the user has committed (clean tree), a full stage-G run into
   `model_test_results/` at builder-default physics.
