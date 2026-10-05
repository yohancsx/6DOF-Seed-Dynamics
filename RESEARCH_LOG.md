# Research log

What was tried, in what order, what it showed, and why the direction changed. The
[README](README.md) says what the code *is*; this file records how it got there.

**How to read an entry.** Each phase gives the question, what was built or changed, the
measured result, and the decision. Commits are the tree the work landed in. Snapshot names
(`model_test_results/<stamp>_<hash>`) identify baseline runs. The snapshot data is archived
outside the repo, but each name carries the git hash that regenerates it.

**Superseded claims are kept, not deleted**, and marked with what replaced them. Several
conclusions below were later corrected, and the corrections are part of the record.

**Reconstructed material is flagged.** Entries before September 2026 are reconstructed
conservatively from commit dates, file-introduction dates and the configuration files of the
time. Where the reasoning was not written down, the entry says so rather than guessing.

---

## Timeline

| date | commit | phase | one line |
|---|---|---|---|
| 2026-06-14 | `de34282`, `104347e` | 0 | Repository and initial seed code |
| 2026-07-13 | `c8dab99` | 0 | Large refactor into the `physics/` layout |
| 2026-07-15 | `d40e371`, `cbe1608` | 0 | First working 6-DOF run; derivations; whole-seed spanwise force added |
| 2026-07-16–21 | `b24056b`–`aae379a` | 0 | Lift/drag multipliers, case-based test suite, first flight-mode tests + classifier |
| 2026-07-23–24 | `7bdb171`, `8862ecc` | 0 | Span torque, CoP migration, attenuation, `Tx`; "FIXED THE DYNAMICS" — all modes appear |
| 2026-08-03–14 | `ba88e49`–`a51df23` | 0 | Mode ablation, mode grid, classifier tuning, grid picker |
| 2026-08-15–26 | `39bdb35`, `6291a28` | 0 | Spanwise twist; moving CoM with twist |
| 2026-09-17 | `af98e91` | 1 | Dual model: frozen `planar` beside new `shape3d` (twist + curvature) |
| 2026-09-20 | `42a6e59` | 1–3 | Physics audit; invented terms cut from `shape3d`; benchmark gate |
| 2026-09-21 | `331ff04` | 4 | Edge drag, added-mass rate, LEV implemented |
| 2026-09-22–23 | `46ab3f1` | 5 | Hou et al. comparison; Kirchhoff pair on by default; 0% map agreement |
| 2026-09-24–26 | `46ab3f1` | 6 | 3D test suite; roll/slide degeneration diagnosed; low-AR plan |
| 2026-09-27 → 10-03 | (uncommitted) | 7 | By-eye review set; paper-metric sign error fixed; six-part classifier |

---

## Phase 0 — Building the 2D→3D strip model (June–August 2026) · *reconstructed*

**Starting point.** The 2D Andersen–Pesavento–Wang (APW) falling-plate model
[[2](README.md#references), [3](README.md#references)], as implemented by Olivia Pomerenk
(`Olivia Code/minimal_imp_Commented.m`), lifted to 6-DOF by slicing the wing into spanwise
strips. Each strip is treated as a 2D APW section, and the strip loads are summed into a
13-state rigid-body model with a discrete nut mass.

**What is known from the record.**

- *07-15.* First working run. The derivation `.tex` files were added the same day. The
  whole-seed **span force** (`computeSpanForce`) arrived in `cbe1608`, the same day. It is
  APW applied a second time in the span–normal plane, because the chordwise strips produce no
  spanwise force by construction. `spanSpinDamping` (`Tx`) and `normalSpinDamping` (`Ty`)
  were in the first working commit.
- *07-16 to 07-21.* Per-strip lift/drag multipliers (`liftMult`, `dragMult`), the case-based
  sweep suite (`runSeedTestSuite`), the first flight-mode suite, and a first classifier.
- *07-23.* One commit ("Added a TON of things") introduced the **span torque**, its
  **migrating span CoP**, the **reduced-frequency attenuation**, the `enableTxDamping`
  switch, `runSpanForceComparison` and `runComMovementTest`.
- *07-24, `8862ecc`.* Commit message: "Dynamics are now able to show all modes! See the
  working inputs." The tuned configuration is recorded in
  `testing/planar/inputs/Working Dynamics Inputs 7-24-26.txt`.
- *08-02.* A second inputs file records the configuration that became the planar default,
  "FULL-minus-geomVelocity": span force, span torque and CoP migration on, `Tx` on,
  attenuation off, `C_span = 0.2`, `C_span_torque = 0.7`. The builder comment says
  `enableSpanGeomVelocity` was dropped because it "had no effect on any mode".
- *08-03 to 08-14.* `runModeAblation` (goal, from its header: "find the SIMPLEST physics
  config ... that still produces all modes"), then the chord × span mode grid, classifier
  calibration against the seven mode-eliciting inputs, and the click-to-animate grid picker.
  Later tuning added a `tiltChaos` gate, a 12 s settling time and a helix-radius fix. That cut
  phase-map "islands" from ~51 to ~31; the rest are real spin-onset bistability.
- *08-15 to 08-26.* Spanwise **twist**: a centred-CoM twisted seed spins up on its own, which
  is shape-driven autorotation with no mass offset. After that, moving-CoM runs with twist.

**Not recorded.** No written record survives of:

- why a whole-seed span force was chosen over other spanwise treatments;
- how `C_span = 0.2` and `C_span_torque = 0.7` were arrived at (beyond "tuned until all
  modes appear");
- what exactly "fixed the dynamics" on 07-24 (the diff touches `computeAeroCoeffs` and
  `seed6DOFODE`).

A **"combined" two-plane model** (a full second, spanwise strip discretization) was
considered and set aside at some point in this period; the reason is not recorded.

**Where it ended.** A model that reproduced every biological mode (7/7 reference inputs),
but only with four empirical whole-seed terms (span force, span torque with migrating CoP,
`Tx`, `Ty`) carrying tuned constants. This is the model frozen as `planar`.

---

## Phase 1 — Dual model and curvature benchmark (2026-09-17 → 09-20)

**Question.** Can real geometry (twist, curvature) supply the out-of-plane forces and moments
that the invented whole-seed terms were standing in for?

**Built.**

- `physics3d/`, a separate non-planar model: `setupSeedShape3D` + `seed6DOFODE3D`, per-strip
  local frames, and curvature-correct mass and inertia.
- The frozen `planar` model beside it. Both call the shared rigid-body core
  `physics/rigidBody6DOF.m` and are selected by `cfg.shapeModel` through `seedRHS`.
- `runCurvatureSuite` / `runCurvatureTest` (symmetric bowl + single side, 0–35° tip
  dihedral, with a separate flat reference row).

**Result.**

- A **symmetric bowl produces no spin at all** (0.00 rad/s throughout), but collapses the
  cone 51°→0.1° and turns the seed into a parachute.
- **Asymmetric** curvature does spin it, up to ~50 rad/s.
- The out-of-plane CoM shift is curvature's signature: +1.40 mm at 35°, matching the hand
  calculation in `testCurvedMassProperties` exactly.

**Decision.** Twist and curvature are complementary mechanisms (spin vs cone/descent). With
geometry available, audit the invented terms.

## Phase 2 — The physics audit: cut the invented terms from `shape3d` (2026-09-20)

**The audit (Sep 2026).** Every invented term existed to compensate for one structural
limitation: a flat, single-plane strip decomposition cannot produce out-of-plane force or
moment. Measured findings:

- **`Tx` double-counts the strips' own roll moment.** A pure-roll state (ω_x only, v = 0,
  CoM centred) isolates the strips' roll moment, which converges on `Tx` as strips refine:
  ratio `1.0141` at 12 strips → `1.0009` at 48 → `1.0001` at 192, holding at `1.0002` across
  ω_x = 5→60. Both reduce to `-ρ·c·CD₂·ω|ω|·S⁴/64`, so roll damping ran at **2× reality**.
- **The span torque is the strips' own normal force re-applied at a fictitious span CoP.**
  `seed6DOFODE` discards `F_span_full(2)` from the force sum as a double count, then crosses
  it with a ~15 mm arm anyway. Over a settled autorotation:
  - `tau_span` is **25.1%** of `|tau_body|`.
  - It splits `6.861e-05` from the discarded force against `1.070e-07` from the force actually
    applied, so ~100% comes from a force that is never applied.
  - The applied span force is `1.641e-05` N, **1.0% of weight**.
  
  This also falsified the code comment claiming strip theory cannot produce roll from
  span-offset loading.
  > *Reinterpreted in phase 6.* At small spanwise incidence this moment turns out to match
  > slender-wing theory to 82–93%. It was mis-derived and mis-shaped, but not arbitrary; see
  > phase 6.
- **The added-mass rate `-Ȧv` is not negligible.** `translationDynamics`' docstring called it
  "tiny for a seed in air". It was finite-difference-verified to `6.2e-09` and runs
  `1.389e-03` N in settled autorotation: 31% of the net force, 82% of weight. Added mass is
  6.3% of `M`, but the rate term is amplified by ω (~19 rad/s).
- **`Ty` is dimensionally a force** (`ρ·R⁴·ω²` is N). That makes it ~20× an honest
  skin-friction estimate, but only 0.1% of the torque budget.

**Changed (in `physics3d/` only; planar frozen byte-identical).**

- Removed from `seed6DOFODE3D`: `Tx`, `Ty`, the span torque, span-CoP migration and the
  attenuation, along with the constants they read (`C_span_torque`, `k0_spanTorque`, `C_Tx`,
  `C_fy`).
- The span **force** was deliberately kept for now, with its honest `cross(r, F)` moment
  (`1.07e-07` N·m).
- `setupSeedShape3D` strips the dead switches, and `buildSeedParams` warns on stale
  overrides.
- Planar gained `enableNormalSpinDamping` (default `true`) so it can be configured down to
  the bare core.
- `testShape3DFlatEquivalence` was re-scoped: bit-identical over the shared core, and
  required to *differ* at full defaults.

**Result: the twist sweep.** Anti-symmetric family, pre-cut → post-cut:

| tip twist | 0° | 5° | 10° | 15° | 20° | 25° | 30° |
|---|---|---|---|---|---|---|---|
| spin, pre (rad/s) | 0.0 | 40.1 | 56.2 | 61.0 | 64.3 | 67.5 | **71.2** |
| spin, post (rad/s) | 0.0 | 50.7 | 81.8 | **86.9** | 65.8 | 60.3 | 54.6 |
| descent, pre (m/s) | 1.66 | 2.00 | 2.76 | 2.83 | 2.87 | 2.93 | 3.02 |
| descent, post (m/s) | 1.66 | 13.66 | 13.42 | 13.47 | 9.24 | 6.91 | 6.09 |
| cone, pre (°) | 51 | 47 | 66 | 68 | 68 | 68 | 68 |
| cone, post (°) | 51 | 88 | 87 | 87 | 78 | 72 | 73 |

**Reading.**

- The twist mechanism does not depend on the cut terms: shape-driven autorotation was never
  an artifact.
- But descent rises 4–5× as the seed rolls edge-on (cone 68°→87°), in converged attractors.
  The phantom moment plus the doubled roll damping had been holding the cone shallow.

**Decision.** The 3D model is not usable quantitatively until real physics replaces what was
cut. That is the argument for phase 4, not against the cuts, which removed provably wrong
terms.

## Phase 3 — Benchmark gate: what survives? (2026-09-20)

Snapshot `2026-09-20_134952_af98e91`. It added a `shape3d/` stage (four twist/curvature
families × 0–35°, centred nut) and ran the grid, CoM and test-suite stages under both models
on the same flat seed.

- **Planar is frozen, verified three ways** against the `a51df23` baseline: byte-identical
  400-cell grid, byte-identical test-suite summary, identical CoM dwell-mode sequences.
- **The cuts mostly moved the model toward real samaras:**

  | | planar | shape3d | delta |
  |---|---|---|---|
  | diving cells | 88 | 56 | **−32** |
  | autorotation cells | 163 | 178 | +15 |
  | spiral / tightSpiral | 83 / 55 | 95 / 62 | +12 / +7 |
  | mean descent (unchanged cells) | 5.27 m/s | 4.20 m/s | **−1.07** |
  | mean cone | 66.3° | 61.2° | **−5.1°** |

  27.5% of cells changed mode, dominated by `diving → autorotation` (×33) and
  `diving → spiral` (×22).
- **The one regression.** Of the seven mode-eliciting inputs, planar scores 7/7 and shape3d
  6/7. The failure is `flutter+spiral`, a nut 0.5 mm off centre (essentially a centred CoM),
  which collapses edge-on to 13.65 m/s at 89.2°.
- **The failure is attitude, not lift magnitude.** The collapsed cases reach genuine terminal
  equilibrium (vertical aero force / weight = `1.00`). The strip angle-of-attack
  distribution is nearly identical across every case: attached ~10%, 20–60° band ~22%,
  bluff ~57%.
- **Diagnosis.** An offset-CoM seed gets its restoring moment from the real weight × arm
  couple, and the cuts helped it. A centred or shape-driven seed has no such couple, and the
  phantom span torque had been standing in for the missing one.

> *Superseded.* "The span force is measurably inert" (≤ `7e-4` relative on every metric) was
> scoped wrong. That seed never slides along its span. In the collapse cases the span force
> was the *only* spanwise resistance: the collapsed seed falls **span-first**, the strips see
> zero in-plane flow, and terminal velocity is set by the spanwise drag term alone.
> `sqrt(2W/(ρ·C_d·A))` with the span force's `CD0·C_span·S·c̄` predicts **13.57** m/s, and
> 13.65 was measured. With no spanwise term the case approaches free fall (82 m/s).

## Phase 4 — Add real physics back (2026-09-21)

Phase 3 reordered the list: the failure is the moment balance, so moment-acting items came
first.

1. **Induced inflow + Prandtl tip loss: dropped before implementing.**
   - *Induced inflow.* Both regimes sit at `V_d/v_h ≈ 23`, far past the windmill-brake
     threshold (|V_d/v_h| > 2), where momentum theory gives `v_i` = 0.19% of descent.
   - *Tip loss.* A centred-CoM wing straddles the rotation axis, so Prandtl's `F` is
     symmetric and adds exactly zero rolling moment to the failing regime. At these inflow
     angles it gives `F = 0.23–0.69` across the whole span, which would cut lift 30–77%.
   - Revisit only if descent brings the rotor back toward `V_d/v_h ≈ 4`.
   - *Note:* this is the rotor tip-loss factor. The planform (Helmbold) lift-slope
     correction is a different thing, and stays open (phase 6).
2. **Finding: the collapse is a chordwise (pitch) problem.** Holding the spanwise offset in
   the collapse band and sweeping the nut chordwise, mid-chord is the *worst* CoM position:

   | nut x/c | CoM, % chord | mode | descent | cone |
   |---|---|---|---|---|
   | −0.25 | 39.1 | fluttering | **2.18** | **49.4** |
   | 0.00 | 50.0 | autorotation | **13.66** | **86.6** |
   | +0.25 | 60.9 | fluttering | **2.19** | **49.4** |

   The best position (39%) is adjacent to the 27–35% of chord behind the leading edge that
   Norberg (1973) reports for flat-plate pitch stability. The governing quantity is the
   chordwise CoP-to-CoM relation, i.e. APW's `l_cp(α)`.
3. **LEV vortex lift: implemented, off by default** (`enableLEV`; `computeLEVForce`,
   `levPlanformConstants`).
   - **Form.** Rezgui et al. (2020), eqs. 3–5, adapting Polhamus:
     `C_L,v = K_v·sin²α·cosα`, `K_v = K_p − K_p²·K_i`.
   - **Constants.** `K_p`, `K_i` come from Helmbold + elliptic induced drag, from the same
     model. Mixing APW's `CL1 = 5.2` with a 3D `K_i` spuriously drives `K_v → 0.04`. Results:
     `K_v` = 1.29 / 2.35 / 2.85 for AR 1.67 / 3.33 / 4.38.
   - **Three findings demoted it as a collapse fix:**
     - Rezgui validate only α = 0–25°, while the term peaks at 54.7°.
     - Snyder & Lamar (1972) put vortex-lift loading at about the potential-flow CoP, so it
       barely touches the pitch balance. The unsourced "forward" option moves the CoP at most
       ~0.035 c, against the ~0.25 c that rescued the collapse.
     - It cannot act span-first, where there is no in-plane flow.
   - **Gate (ours, not Rezgui's).** `G = 1/(1+(Ro/Ro_crit)^p)`. Lentink & Dickinson (2009)
     publish no critical Rossby number, so `Ro_crit = 3`, `p = 4` are anchors.
   - **Rossby definition.** The first gate read `Ro = r/c` geometrically. That gave a
     parachuting seed that never spins `G ≈ 1`, breaking reference modes (5/8 preserved).
     Replaced by the **kinematic** `Ro = |v_ip|/(Ω·c)`, `Ω = |ω × ŝ|`, Rossby's general
     ratio. It reduces exactly to `r/c` under pure revolution (verified ~1e-16) and preserves
     8/8 reference modes; the geometric form is kept as a switch.

     | | LEV off | `'geometric'` | `'kinematic'` |
     |---|---|---|---|
     | `parachute` mode | ✓ | → gliding | ✓ |
     | `autorotation` mode | ✓ | → spiral | ✓ |
     | `flutter+spiral` (collapse) mode | fluttering | → autorotation | fluttering |
     | tight-spiral descent | 1.44 | 1.23 | 1.30 |
     | twist 10° vertical spin | 60.1 | 98.0 | 60.0 |
     | **reference modes preserved (of 8)** | — | **5** | **8** |

   - **Circularity.** The gate is half open at `Ω* = |v|/(Ro_crit·c)`, ~89 rad/s at 4 m/s for
     a 15 mm chord, while the model autorotates at 17.8 rad/s. The lift deficit makes the
     seed fall fast, and falling fast denies it the lift.
   - *(An earlier worry that the model cannot tell leading from trailing edge was wrong: the
     leading edge is `sign(v_c)`, which `computeAeroCoeffs` already branches on.)*
4. **Added-mass rate `-Ȧv`: implemented** in the shared core behind `enableAddedMassRate`.
   - Initially **off by default** while its literature support was reviewed.
   - Check A5 shows that with it on, the translational equation *is* Kirchhoff's body-frame
     form, the `(m+m2)·ω·vyp` / `−(m+m1)·ω·vxp` terms of `minimal_imp` (lines 208–209).
     With it off, the equation is off by up to 140% on the same states.
   - *Turned on by default in phase 5.*
5. **Nut form drag.** Judged least relevant (a centred nut acts near the CoM). Still open.
6. **Edge crossflow drag: implemented, on by default** (`computeEdgeDrag`,
   `F = −½ρ·C_d·(t·c̄)·|v_s|·v_s·ŝ`, `C_d = 1.2`). It **replaced `computeSpanForce`**, which
   was retired from `shape3d` along with `enableSpanForce` and `enableSpanGeomVelocity`.
   > *My prediction that it would be inert was wrong.* In the collapse cases it is the only
   > spanwise resistance and sets terminal velocity exactly: predicted 8.76 m/s, measured
   > 8.76 m/s.

**Benchmark after items 4 and 6.** Snapshot `2026-09-21_103642_6940ceb`, run with the rate
term on.

- Planar is still byte-identical to `a51df23`.
- **Against the previous 3D grid:** 43/400 cells changed. Cells descending faster than 10 m/s
  went from 33 to 0, now bounded by honest tip drag (~8.8 m/s) instead of the span force
  (~13.6 m/s).
- **Against planar:** diving −26, autorotation +25, and mean descent 5.24 → 3.88 m/s.
- Reference inputs still 6/7: `flutter+spiral` collapses to `diving` at 8.76 m/s.

## Phase 5 — Match the Hou et al. (2025) phase map (2026-09-22 → 09-23)

**Target.** Fig. 2a of Hou et al. [4]: an *experimental* four-mode map over CoM position,
with these modes:

- autorotation (AR)
- spiral tumbling (ST), continuous or segmented
- chaotic (CH)
- falling (FA)

Their plate is 60 × 20 mm (AR 3.0, 92.8 mg), carrying a 102.9 mg weight on the long axis and
a 51.4 mg weight on the short axis. The goal was to reproduce the map's topology and
descent-speed ordering with minimal added physics, not to fit it.

**Decoding their axes.** The caps are structural, which is what confirms the reading:

| theirs | is | cap | ours |
|---|---|---|---|
| `x_c = m_h·x′/m_tot`, `a = L/2` | CoM offset along the **long** axis (our span) | 0.416 | `2µ·spanFrac` |
| `y_c = m_l·y′/m_tot`, `b = W/2` | CoM offset along the **short** axis (our chord) | 0.208 | `2µ·chordFrac` |

Our existing grid ran both axes to 1.30, so their whole map was the bottom-left ~5% of ours.
Only 2 of our 7 reference inputs fall inside it.

**Tier 0: measurement, no new physics.** Built:

- `paperSeedConfig`: their plate, masses, release (long axis horizontal, from rest), and a
  thickness derived from their areal density.
- `runPaperModeGrid` + `reportPaperModeGrid`: a separate grid, with reporting split out so a
  reporting bug cannot lose the integrations.
- `computePaperMetrics` + `classifyPaperMode`: their vocabulary, reported beside ours.
- `runPaperAnchorCases`: the five published configurations with targets.
- `paperModeReference`: a hand-digitised Fig. 2a.

**Classifier defects found by `testPaperMetrics`** (synthetic motion with known answers,
8 cases):

1. **`ω_z` read directly.** A pure tilted revolution at 40 rad/s showed 7.9 rad/s of phantom
   tumbling and was labelled CST. Fixed with self-rotation `ψ' = ω_z − φ'·sinθ`.
2. **Segmented ST labelled AR**, because direction reversals zero the net rotation. Fixed with
   the largest one-way rotation between reversals.
3. **CH required non-convergence**, so a ±44° attitude swing with steady descent read as FA.
   Fixed to key on orientation instability.

A flutter label `FL`, deliberately outside their vocabulary, was added after animations
showed rocking at ψ' = 45 rad/s with 218 reversals and 0.2 turns.

> *Superseded.* An earlier reading of the grid reported 28% agreement with ST absent. It came
> from the classifier with defect 1.

**Tier 1 screen** (R = added-mass rate, M = added-mass moment, L = LEV; five anchors):

| combo | AR θ | CST θ | RMS `V_d/scale − 1` | FA fastest |
|---|---|---|---|---|
| (baseline) | −18.3 ± 28.0 | −56.1 ± 0.2 | 0.89 | no |
| `R` | −19.7 ± 26.1 | −58.7 ± 0.0 | 2.10 | no |
| `L` | −27.7 ± 6.1 | −73.8 ± 18.0 | 4.85 | no |
| `M` | −32.6 ± 0.9 | −0.1 ± 0.1 | 0.24 | **yes** |
| **`RM`** | **−29.6 ± 1.0** | −0.1 ± 0.1 | **0.20** | **yes** |
| `RML` | −27.9 ± 3.7 | −0.1 ± 0.1 | 0.24 | **yes** |
| *paper* | *−11.4 ± 4.2* | *−38.2 ± 2.3* | *~0* | *yes* |

**How the Munk moment was found.** An earlier Tier 1 idea, "apply edge drag at the windward
tip", was *wrong*: for a flat seed the edge force is parallel to its own moment arm, so
relocating it cannot change `r × F`. Chasing that found the real gap. Kirchhoff's equations
are a **pair** (Lamb, ch. VI), and only the translational half had been implemented. APW
eq. (6.3) carries the rotational half explicitly as `(m₁₁ − m₂₂)v_x′v_y′`. Implemented as
`enableAddedMassMoment` (`−v × (A·v)`); checks M1–M6 verify the 2D limit against APW (6.3)
to `2.2e-16`.

**Decision: R + M on by default in `shape3d`.** The rule is that an added-mass term belongs in
by default when it is first-principles *and* present in the reference model (APW or the 2D
baseline). Both halves of Kirchhoff's pair are. They moved descent from 1.2–2.4 to
1.0–1.4 × `sqrt(σg/ρ)` and collapsed the AR case's cone wander from ±28.0° to ±0.9°. Planar
keeps both off, to stay byte-identical. LEV stays off: near-neutral with M on, harmful alone.

> *Superseded numbers from this phase.*
> - The pre-Tier-1 `V_d/sqrt(σg/ρ)` of 3.0–6.5 (our seed) became 1.04–1.36 on the paper
>   plate with R+M.
> - "Their FA mode is our collapse" did not survive the grid: in their FA region we produce
>   ST.
> - "`I_long` is ~5× too low" was computed on our 15 mm seed. On their plate the lumped nut
>   is ~1.6× low at the top of the `y_c/b` range, too small to be the cause.

**Measured grid.** Snapshot `2026-09-23_184153_331ff04`, 13×8 over their window, shape3d
defaults. **Agreement with digitised Fig. 2a: 0% (0/104).**

```
OURS                            THEIRS (digitised)
y=0.200 | LLLLscccccccc         y=0.200 | XAAAAAAAAAAAA
y=0.146 | LLLLLsccccccc         y=0.146 | XXXXXAAAAAAAA
y=0.091 | LLLLLLscccccc         y=0.091 | TTTXXXXXXFFFF
y=0.037 | LLLLLLLcccccc         y=0.037 | TTTTTTTXXFFFF
        L=flutter  c=CST  s=SST         T=ST  X=CH  A=AR  F=FA
```

| they say | cells | we say |
|---|---|---|
| ST | 25 | flutter ×25 |
| CH | 27 | flutter ×17, ST ×10 |
| AR | 36 | **ST ×32**, flutter ×4 |
| FA | 16 | **ST ×16** |

**Reading.** Descent speed is right everywhere (1.04–1.36 × `sqrt(σg/ρ)`). Every failure is
an attitude failure: **the model cannot hold a fixed attitude about its long axis.** It
either flutters, or tumbles while revolving at a steady −20 to −30° cone. Both missing modes
are defined by the *absence* of that flipping, so one deficit removes both. Their AR is
`θ = −11.4 ± 4.2°` with no tumbling; ours is `−29.6 ± 1.0°` at 0.97 tumbles per revolution.

**Ruled out by measurement:**

- a mis-fitted revolution rate (ground-track and span-azimuth fits agree to four significant
  figures);
  > *Superseded (phase 7).* This check could not see the error that mattered. Both fits
  > measured angle about −Y, while the body angular velocity is about +Y, so they agreed
  > with each other and were both opposite to the body rotation. That sign error made every
  > steady revolver read as spiral tumbling. "They say AR → we say ST ×32", and this
  > phase's 0% agreement, need re-measuring with the fix.
- the `I_long` gap (1.6×, see above);
- release attitude and strip count (checked on the anchors).

**Open question left.** Why does the chordwise CoM offset produce so little pitch stiffness?
Their AR needs `y_c/b ≳ 0.12`; our map barely responds to `y_c/b`.

**Known ceiling.** Their ST lift peak and AR mode involve tip-vortex structures with memory
(Ω-shaped tubes, rib-like vortices) that quasi-steady strip theory cannot carry.

## Phase 6 — 3D test suite and the roll/slide diagnosis (2026-09-24 → 09-26)

**Built.**

- `runSeedTestSuite3D`: the one-click 3D suite. It runs twist and curvature (with
  animations), the nine 2D-parity sweep groups, the moving-CoM scenarios, our-seed mode grid,
  and the paper comparison.
- `Seed_Dynamics_ODE_Test_3D`, the single-drop harness for the 3D model.
- `animateCaseVideo` and `buildComPath` as shared helpers.
- A fix to `pickModeGridRuns`, which hardcoded the planar RHS and would have animated
  different physics from the one that labelled a 3D cell.
- The twist sweep now releases level, so the spin-up is the shape's and not the release
  tilt's.

**Observed** (hand-checking suite animations):

1. Cells are misclassified, partly through threshold tuning and partly because the model
   produces **modes real seeds do not**. The clearest is an edge-on mode with the span
   pointing straight down.
2. Without the planar patches, the seed **degenerates by rolling and sliding spanwise**.

**Diagnosis.** Both involve flow with a spanwise component, which strip theory discards by
construction. That is exact for an infinite wing; at AR ≈ 3, tip effects are first-order.

- **Roll + slide.** The strips' rolling moment does not depend on sideslip at all
  (C_lβ = 0). A rolled plate slides under its tilted normal force and nothing rolls it back.
  A real thin plate also has little spanwise drag. What stops it is sideslip loading the
  windward tip, which banks the plate. So the missing piece is a **moment, not more drag**.
- **Span-down mode.** Near end-on, the strip force is uniform along the span, so there is no
  moment to leave. A real slender plate has lift at its leading tip, linear in the tilt,
  which makes end-on unstable.

**Finding: the planar span patch was mostly real physics.** Its torque for near-spanwise flow,
against Jones's slender-wing moment (lift `(π/2)·q·c²·β` at the leading tip), on the
50 × 15 mm seed:

| spanwise incidence β | 1° | 2° | 5° | 10° | 20° |
|---|---|---|---|---|---|
| patch / Jones | 0.93 | 0.92 | 0.82 | 0.52 | 0.17 |

It falls off past ~10° because it borrowed APW's 2D stall (14°), which a slender plate does
not have. So the patch had the right small-angle physics, a borrowed derivation, and the
wrong large-angle shape. Second instance: Helmbold's attached lift slope is 3.4–3.6/rad at
AR 3–3.3, against APW's `CL1 = 5.2`, so the strips over-predict attached lift on these
plates by ~50%.

**Decisions** (2026-09-26):

- **Rule for new terms.** A term's **form** comes from theory and its **magnitude** from
  canonical *static* published data (fixed plates, tunnel or DNS). That is the same standing
  the APW coefficients have. Free-flight data (Hou et al., our own trials) is validation only,
  never tuning.
- **Order of work.**
  1. Fix the classifier (non-physical labels; a hand-labelled reference set with holdout).
  2. Build a virtual wind tunnel comparing static forces and moments against published
     low-AR plate data.
  3. Add low-AR physics one piece at a time, each accepted only if it improves the static
     comparison: sideslip roll moment, slender-wing lift, finite-AR attached lift, then
     side-edge vortex lift.
  4. Validate against trials with pre-registered predictions and mode probabilities.
- The open items are in the README roadmap.

## Phase 7 — Fixing the classifier (2026-09-27 → 10-03)

**Question.** Is the classifier labelling the motion correctly? Every comparison, including
the paper map, runs through it.

**Built.**

- `buildClassifierReviewSet`: a 6×6 nut grid on our seed under both models, 72 runs.
  Shuffled IDs; blind animations with no mode colours or labels; static summaries with
  cumulative rotation; saved trajectories; and a labelling spreadsheet.
- `MODE_DEFINITIONS.md`, the definitions file.

**By-eye labels** (review set `2026-09-27_165346_3727986`). These were the observations:

- 12 different label strings, plus 13 rows with notes only.
- The labels were already two-part: rotation about the span (tumble, flutter) plus what the
  path does (dive, glide, spiral). The inconsistency came from naming both in one word:
  - one coned-revolution cluster got 4 different names;
  - one flutter-while-revolving cluster got 4 names.
- "Dive" covered three different states:
  - chord edge first, span level (6.1 m/s);
  - edge-on with the span tilted 53–64° (9–11 m/s);
  - span-first, span 82–85° below horizontal (11–14 m/s). This is the non-physical end-on
    mode, which on video just reads as a dive.
- A note on R003 gave a usable rule: autorotation when the yaw axis is within about a span
  of the seed, otherwise a spiralling glide. Helix radius over span splits the runs with an
  empty gap between 0.25 and 1.5.

**Bug found by the labels: sign error in the paper metrics.**

- *Symptom.* On every steady revolver, the self-rotation rate was exactly twice the raw span
  rate (R035: 91.7 vs 45.9 rad/s) instead of zero.
- *Cause.* `computePaperMetrics` measured the orbital revolution angle as `atan2(+dZ, dX)`,
  the angle about −Y. The body angular velocity is about +Y. So subtracting the revolution's
  projection doubled it instead of cancelling it.
- *Why the test missed it.* `testPaperMetrics` built its synthetic orbit as
  `[R cos φ, ·, +R sin φ]` while turning the body by `Ry(φ)`. That flip matched the
  metric's, so it passed.
- *Fix.* The metric now uses `atan2(−dZ, dX)`, and the synthetic orbit is generated by the
  same rotation matrix as the attitude.
- *Effect.* R035 went from 91.7 rad/s of "tumbling" to 0.00, and the paper classifier now
  calls it AR. **Every paper-classifier result before 2026-10-03 is affected**, including the
  phase-5 grid.

**Decisions** (2026-10-03):

- Autorotation vs spiral glide at a revolution radius of **1 span**.
- End-on (non-physical) at a span angle of **70°** from horizontal.
- Glide vs parachute at a glide ratio of **0.25**.
- **Wobbling spin and chaotic are non-physical.**
- Tumbling without revolution is `tumble`. With revolution it is `tightSpiralTumble` when
  the radius is under a span, and `spiralTumble` otherwise.

**Built: the six-part classifier** (`classifySeedMode`, `defaultSeedModeThresholds`,
`computeSeedModeMetrics`). It judges six things independently and builds the name from
them, so the same motion always gets the same name:

- flip about the span;
- revolution about the vertical;
- span attitude;
- path;
- regularity;
- physicality.

Two measurement choices remove the failure modes found above:

- **Flip is measured from attitude:** the plate normal's angle about the span, referenced
  to up. It never touches the revolution rate.
- **Revolution is measured from the span axis's heading**, not the body spin about the
  vertical. A tumble about a tilted span projects onto the body spin, but leaves the
  heading fixed.

Two first-draft rules were corrected on the review data:

- An edge slide needs a *level* span. A steep span falling straight down also lines up with
  the velocity.
- Chaotic also triggers when the span-angle spread exceeds 25°. R030, at 31°, had been
  called a wobbling spin.

**Result.**

- **Synthetic regression:** 17/17 cases pass (8 paper-metric, 9 new-classifier, including
  the tilted-tumbler trap).
- **On the 72 review runs:** the labels match the distilled scheme. The remaining
  disagreements with the eye labels are listed in `MODE_DEFINITIONS.md` section 3 for
  re-watching.

---

## Decided against or superseded

| item | status | why | where |
|---|---|---|---|
| Tune the span torque (`C_span_torque`) | decided against | Phantom moment in its lumped form; to be replaced by a theory-shaped term (phase 6) | phases 2, 6 |
| Reduced-frequency span-torque attenuation | cut from `shape3d` | Existed to suppress the lumped span CoP's spin-frequency oscillation; a per-strip formulation should not need it (to be checked) | phases 2, 6 |
| Unsteady CoP lag (first-order state) | superseded | Was the principled replacement for the attenuation, which is gone from `shape3d` | phase 2 |
| Crossflow-drag-only span-force variant | superseded | Edge crossflow drag implemented instead | phase 4 |
| Keep the span force in `shape3d` | retired | Replaced by edge drag; its drag level had been setting collapse terminal velocity | phases 3, 4 |
| Induced inflow + Prandtl rotor tip loss | dropped | `v_i` ≈ 0.19% of descent; symmetric tip loss adds no roll moment | phase 4 |
| LEV as the collapse fix | demoted | CoP effect ~7× too small; cannot act span-first; kept as a switch for lift magnitude | phase 4 |
| Edge drag relocated to the windward tip | wrong | Force parallel to its arm; `r × F` unchanged. Led to the Munk moment | phase 5 |
| Prune dead aero constants when planar retires | not planned | Planar is kept frozen as the reference; the constants stay | phase 2 |
| Fix the `(π−|α|)` branch quirks | not worth it | Cost ≤ 0.7% of peak CT, 0.2% of peak CD, measured | phase 0 port |
| Full "combined" two-plane spanwise re-discretization | set aside (reason unrecorded) | The phase-6 low-AR plan revisits the same spanwise-flow problem | phases 0, 6 |
