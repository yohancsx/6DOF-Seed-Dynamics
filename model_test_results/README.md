# Model test results

A record of **what the model does** at a given state of the code and physics. When the
model changes, run the suite again and compare the new run folder against an earlier one.

## Suite runs

Every run of [`testing/shape3d/runSeedTestSuite3D.m`](../testing/shape3d/runSeedTestSuite3D.m)
writes one folder:

```
model_test_results/
  <YYYY-MM-DD_HHMMSS>_<githash>[-dirty]_3D_<label>_<profile>/
    summary_3D.txt              everything, opening with the regression-test verdict
    regression_3D.txt / .csv    stage 0: every pass/fail check
    tracking_3D.txt             this run's tracked numbers
    twist_3D/  curvature_3D/  com_movement_3D/      stages A, B, D (+ videos, full profile)
    <group>/                    stage C, the nine sweep groups (.fig/.png/.mat, overlays)
    mode_grid_3D.mat / .png     stage E, our seed's mode grid
    mode_grid_picks/            videos made from that grid with testing/tools/pickModeGridRuns.m
    paper_mode_grid_3D.mat      stage F, the Hou et al. comparison
    paper_comparison_3D.png     ours vs the digitised Fig. 2a, disagreements outlined
    manifest_3D.mat             cfg, switches in force, git, tracked structs
  tracking_history_3D.csv       one row per run (merged, so new columns are safe)
```

Read the history with `readtable(file, 'ReadVariableNames', true, 'Delimiter', ',', 'TextType', 'string')`.
MATLAB's auto-detection can mistake its header for data.

- **Labels:** the folder's `<label>` names the physics configuration (the `sw` overrides).
- **Profiles:** `<profile>` is `full`, for baselines (everything, with videos), or `quick`, a
  fast look after a change (coarse grids, no videos). Compare runs only within one profile.
- **Older names:** runs before 2026-10-15 are named `..._3D_<label>` (no profile) and are full
  runs.
- **Older snapshots:** folders named `<timestamp>_<githash>` with no `_3D_` were made by the
  now-archived `generateModelBaseline`, which ran the planar and 3D models side by side.

## Other folders

| folder | contents |
|---|---|
| `planar_final_<hash>/` | The **final record of the archived planar model**: its seven named modes, its test suite and its mode grid, run at the last commit where it was live (git tag `planar-final`). |
| `fixtures/` | Golden references for `testing/shape3d/testRefactorInvariance.m`: built seeds, RHS samples and trajectories, pinned bit for bit around a refactor. |

## Why it's git-ignored

Runs contain videos and full trajectory `.mat`s and add up fast, so the outputs are **not**
uploaded to GitHub (see this folder's `.gitignore`). Only the folder and this README are
tracked. The **git hash** in each folder name and manifest pins the exact code that produced
it, so any run can be regenerated from source. Older snapshots cited in
[`RESEARCH_LOG.md`](../RESEARCH_LOG.md) are archived outside the repo.
