# Model test results — baseline snapshots

A record of **what the model does** at a given state of the code/physics, kept as
a regression baseline. When the model changes, regenerate a snapshot and diff it
against an earlier one to see what moved.

Each snapshot is a timestamped, git-stamped subfolder produced by
[`testing/baselines/generateModelBaseline.m`](../testing/baselines/generateModelBaseline.m):

```
model_test_results/
  <YYYY-MM-DD_HHMMSS>_<githash>/
    manifest.txt / manifest.mat   # git commit + dirty flag, MATLAB version, seed
                                   # geometry, BOTH models' physics switches, aero +
                                   # classifier settings, every stage's config/inputs
    mode_grid/     mode_grid_3d/   # coarse flight-mode phase map, planar / shape3d,
                                   # same flat seed (modeIdx + metrics + .png + ascii)
    com_movement/  com_movement_3d/  # the 3 moving-CoM scenarios: trajectories,
                                   # per-dwell mode classification, animations
    test_suite/    test_suite_3d/  # the sweep battery (per-case .mat/.png, overlays,
                                   # summary.txt), same as runSeedTestSuite output
    shape3d/                       # twist + curvature families (shape3d, centred nut)
    paper/                         # Hou et al. (2025) comparison: the five anchor
                                   # cases + the mode grid over their window, both
                                   # classifiers per cell, agreement with Fig. 2a
```

Diff a snapshot's planar and shape3d stages with
[`testing/baselines/compareBaselineModels.m`](../testing/baselines/compareBaselineModels.m).

## Why it's git-ignored

Snapshots contain videos and full trajectory `.mat`s and add up fast, so the run
outputs are **not** uploaded to GitHub (see this folder's `.gitignore`). Only the
folder and this README are tracked. The **git hash** in each snapshot's name and
`manifest` pins the exact code (including the aero constants) that produced it, so a
snapshot can be regenerated from source even though the data isn't committed. Older
snapshots cited in [`RESEARCH_LOG.md`](../RESEARCH_LOG.md) are archived outside the repo.

## Regenerating

Run `testing/baselines/generateModelBaseline.m`. It adds the paths it needs and writes a
fresh snapshot folder here.
