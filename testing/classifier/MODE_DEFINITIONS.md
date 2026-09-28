# Flight-mode definitions

The single reference for what every mode label means, in words and mathematically. It
serves three purposes:

- **By-eye labelling.** Use the *eye vocabulary* (section 4) when filling the `label_eye`
  column of a classifier-review spreadsheet
  (`testing/classifier/buildClassifierReviewSet.m`).
- **Record of the current code.** Sections 2–3 transcribe exactly what the two existing
  classifiers do, including their known defects.
- **Growing list.** When you see a new mode, add it to section 4 using the template in
  section 5. Give it a words definition first; the math can come later.

Body axes: `x` = chord, `y` = plate normal, `z` = span. World `Y` is up. `R(t)` is the
body→world rotation.

---

## 1. Shared measurements

Every metric is computed over a **settled window** `W = { t : t ≥ f·t_end }`, with
`f = windowStartFrac = 0.5` by default. For a 12 s drop, that is the last 6 s. Averages are
over the ode45 **samples** in `W`, which are not uniformly spaced in time.

Unit vectors in the world frame:

| symbol | definition | meaning |
|---|---|---|
| $\hat n$ | $R\,\hat e_y$ | plate normal |
| $\hat s$ | $R\,\hat e_z$ | span (long) axis, $+z$ = nut side |
| $\hat Y$ | $(0,1,0)$ | world vertical |

### From `computeTrajectoryMetrics` (used by the current flight-mode classifier)

| metric | definition | units |
|---|---|---|
| `descentSpeed` | $-\overline{v_Y}$ | m/s, + down |
| `glideRatio` | horizontal displacement / height lost, over $W$ | – |
| `helixRadius`, `helixValid` | circle fit to the horizontal ground track $(X,Z)$ | m |
| `coneAngleDeg` $\bar\gamma$ | mean of $\gamma=\arccos\lvert\hat n\cdot\hat Y\rvert$: 0° broadside, 90° edge-on | deg |
| `tiltStd` $\sigma_\gamma$ | std of $\gamma$ over $W$ | deg |
| `verticalSpinMag` $\Omega_V$ | $\overline{\lvert (R\,\boldsymbol\omega)\cdot\hat Y\rvert}$ | rad/s |
| `spanwiseSpin` $\Omega_s$ | $\overline{\lvert\omega_z\rvert}$ (raw body $z$ rate) | rad/s |
| `tumbleFrac` | $\lvert\int\omega_z\,dt\rvert / \int\lvert\omega_z\rvert\,dt$ | 0..1 |
| `converged` | across the two halves of $W$: relative change in descent < 0.20 **and** in $\Omega_V$ < 0.60 (0.1 rad/s floor) | bool |

### From `computePaperMetrics` (used by the paper classifier)

| metric | definition | units |
|---|---|---|
| `spanAxisTiltDeg` $\theta$, `spanAxisTiltStd` $\sigma_\theta$ | mean / std of $\theta=\arcsin(\hat s\cdot\hat Y)$, the long axis from horizontal ($<0$: nut tip down) | deg |
| `revolutionRate` $\dot\phi$ | slope of the unwrapped ground-track angle about the fitted helix centre; NaN if no valid helix | rad/s |
| `selfRotationRate` $\overline{\lvert\dot\psi\rvert}$ | $\dot\psi=\omega_z-\dot\phi\sin\theta$: rotation about the long axis with the revolution's projection removed | rad/s |
| `maxTurnsOneWay` | largest $\lvert\int\dot\psi\,dt\rvert/2\pi$ between direction reversals | turns |
| `tumbleReversals` | reversals of $\dot\psi$: Schmitt trigger at 25% of peak on $\dot\psi$ smoothed over 2% of $W$ | count |

---

## 2. Current vocabulary A — `classifyFlightMode` (our seed)

Thresholds (`defaultModeThresholds`):

| threshold | value |
|---|---|
| `spinLo` | 5 rad/s |
| `vSpinAuto` | 10 rad/s |
| `tiltSteady` | 5° |
| `tiltChaos` | 30° |
| `helixTight` | 0.05 m |
| `glideHi` | 1.0 |
| `coneEdge` | 45° |
| `coneBroad` | 20° |

Let $\text{spinning} = \Omega_V>$ `spinLo` $\lor\ \Omega_s>$ `spinLo`. The rule tree is
applied in order; the first match wins.

| # | label | rule | in words |
|---|---|---|---|
| 1 | `chaotic` | $\lnot$converged $\land\ \Omega_V<10 \land \sigma_\gamma>30°$ | unsettled, no coherent spin, attitude swinging |
| 2 | `gliding` | $\lnot$spinning $\land$ glideRatio $>1$ | no rotation, travels sideways more than down |
| 3 | `diving` | $\lnot$spinning $\land\ \bar\gamma>45°$ | no rotation, edge-on, fast |
| 4 | `parachuting` | $\lnot$spinning $\land\ \bar\gamma<20°$ | no rotation, broadside, slow |
| 5 | `undetermined` | $\lnot$spinning, otherwise | no rotation, intermediate attitude |
| 6 | `autorotation` | $\Omega_V>10 \land \sigma_\gamma<5°$ | spins about the vertical at a steady attitude |
| 7 | `tightSpiral` | $\Omega_V>10$, valid helix with $R<0.05$ m | spins about the vertical, attitude flipping, tight helix |
| 8 | `spiral` | $\Omega_V>10$, otherwise | as above, wide helix |
| 9 | `fluttering` | spinning, $\Omega_V\le 10$ | rotating, but not about the vertical |

Upstream of the tree, runs with non-finite states are labelled `chaotic`, and runs whose
integration errors are labelled `failed`.

**Known defects:**

- **Rule 9 lumps two different motions.** `fluttering` covers both rocking (never a full
  turn) and continuous tumbling (whole turns about the span). `tumbleFrac` is computed but
  unused.
- **`spinning` reads raw $\lvert\omega_z\rvert$.** A seed revolving at a tilt projects its
  revolution onto the span axis. That is the defect `testPaperMetrics` case A caught in the
  paper classifier; this classifier still has it.
- **Attitude is the plate normal $\gamma$, not the span axis.** So rule 3 `diving` cannot
  tell edge-on-chord-down (a real dive) from **span-down** (the non-physical end-on mode):
  both have $\gamma\approx 90°$.
- **Thresholds were calibrated on 7 planar reference inputs only.**

---

## 3. Current vocabulary B — `classifyPaperMode` (Hou et al. 2025)

Thresholds (`defaultPaperThresholds`):

| threshold | value |
|---|---|
| `tumbleTurnsST` | 1.0 turn |
| `tumbleRateST` | 5 rad/s |
| `tiltStdSegmented` | 7° |
| `vSpinAR` | 10 rad/s |
| `tiltStdAR` | 10° |
| `tiltStdChaos` | 25° |

Let:

- $\text{tumbling} = \text{maxTurnsOneWay}>1 \land \overline{\lvert\dot\psi\rvert}>5$
- $\text{revolving} = \Omega_V>10$
- $\text{steady} = \sigma_\theta<10°$

The rule tree is applied in order; the first match wins.

| # | label | rule | in words |
|---|---|---|---|
| 1 | `SST` | tumbling $\land$ (reversals $>0 \lor \sigma_\theta>7°$) | segmented spiral tumbling: whole turns about the long axis, direction reverses |
| 2 | `CST` | tumbling, otherwise | continuous spiral tumbling at a steady cone |
| 3 | `CH` | $\sigma_\theta>25°$ | attitude never settles |
| 4 | `CH` | $\lnot$converged $\land\ \sigma_\theta>7°$ | unsettled with an unsteady cone |
| 5 | `AR` | revolving $\land$ steady | autorotation: revolves about the vertical, steady long-axis tilt, no tumbling |
| 6 | `CH` | revolving, otherwise | rotating but the cone never settles |
| 7 | `FL` | $\overline{\lvert\dot\psi\rvert}>5 \land$ maxTurnsOneWay $<1$ | flutter: fast rocking about the long axis, never a full turn (not a paper mode) |
| 8 | `FA` | otherwise | tilted fall without sustained rotation of either kind |

**Known defects:**

- **ST requires no revolution.** A tumbler that does not revolve (2D-style tumbling) is
  labelled `CST`/`SST`.
- **FA is a catch-all.** It includes the non-physical span-down fall.
- **Validated on synthetic motion only** (`testPaperMetrics`), not against eye labels.

---

## 4. Eye vocabulary — use these in `label_eye`

One label per run, judged on the **settled** behaviour (the video shows the classifier
window). If two apply, pick the dominant one and note the other in `notes_eye`. Paper
equivalents are given where they exist. The math column is a *proposed* definition for the
new classifier; blank means not yet defined.

| label | in words | paper | proposed math (draft) |
|---|---|---|---|
| `gliding` | no sustained rotation; moves sideways at least as much as down, at a shallow attitude | – | no rotation (below); glideRatio $>1$ |
| `diving` | no sustained rotation; falls fast along a steep straight path, **chord** leading edge down | FA (partly) | no rotation; $\bar\gamma>45°$ and $\lvert\theta\rvert<45°$ |
| `parachuting` | no sustained rotation; falls broadside, slowly | – | no rotation; $\bar\gamma<20°$ |
| `falling` | no sustained rotation; steady large tilt, heavier side down, fast | FA | no rotation; $\lvert\theta\rvert$ large and steady |
| `fluttering` | rocks back and forth about the span axis, **never completes a turn** | FL | $\overline{\lvert\dot\psi\rvert}>5$, maxTurnsOneWay $<1$, $\Omega_V$ small |
| `tumbling` | turns end-over-end about the span axis, **whole turns**, no revolution about the vertical | – | maxTurnsOneWay $\ge1$, no valid helix / $\lvert\dot\phi\rvert$ small |
| `spiralTumbling` | tumbles end-over-end **while** the path spirals about a vertical axis | CST / SST | tumbling $\land$ revolving |
| `autorotation` | revolves about a vertical axis at a **steady** attitude, no flipping | AR | revolving, $\sigma_\theta$ small, maxTurnsOneWay $<1$ |
| `chaotic` | no sustained pattern; attitude swings irregularly | CH | $\sigma_\theta$ large, no periodicity |
| `endOn` ⚠ | **non-physical:** falls with the span axis (near) vertical, span-first | – | $\lvert\theta\rvert>70°$ held over $W$ |
| `edgeSlide` ⚠ | **non-physical:** slides along its own span, plate edge-on | – | $\lvert\hat s\cdot\hat v\rvert>0.5$, $\gamma\approx90°$ |
| `unsure` | can't tell; say why in `notes_eye` | – | – |

⚠ = a mode real seeds do not show. The model produces these because strip theory discards
spanwise flow (research log, phase 6). They count as **disagreements** in any score, and
their share of a grid is a tracked metric that should go to zero.

"No rotation", as a draft: $\Omega_V<5$ and $\overline{\lvert\dot\psi\rvert}<5$ rad/s.

---

## 5. Template for a new mode

Copy into the table in section 4:

```
| `newLabel` | what you SEE, in one sentence (attitude, rotation, path) | paper equivalent or – | leave blank, or a first guess |
```

Then add a line below describing one example run (review folder + `run_id`) where it
appears, so the definition can be checked against footage later.

Examples seen:

- *(none yet)*
