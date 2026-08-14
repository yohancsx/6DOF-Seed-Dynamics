# Model test results — baseline snapshots

A record of **what the model does** at a given state of the code/physics, kept as
a regression baseline. When the model changes (e.g. the planar → full-3D work),
regenerate a snapshot and diff it against an earlier one to see what moved.

Each snapshot is a timestamped, git-stamped subfolder produced by
[`testing/generateModelBaseline.m`](../testing/generateModelBaseline.m):

```
model_test_results/
  <YYYY-MM-DD_HHMMSS>_<githash>/
    manifest.txt / manifest.mat   # git commit + dirty flag, MATLAB version, seed
                                   # geometry, physics switches, aero + classifier
                                   # settings, and every stage's config/inputs
    mode_grid/                     # coarse flight-mode phase map (modeIdx + metrics
                                   # + map .png + ascii); re-integrable via saved cfg
    com_movement/                  # the 3 moving-CoM scenarios: trajectories,
                                   # per-dwell mode classification, animations
    test_suite/                    # the sweep battery (per-case .mat/.png, overlays,
                                   # summary.txt) — same as runSeedTestSuite output
```

## Why it's git-ignored

Snapshots contain videos and full trajectory `.mat`s and add up fast, so the run
outputs are **not** uploaded to GitHub (see this folder's `.gitignore`). Only the
folder + this README are tracked. The **git hash** recorded in each `manifest`
pins the exact code (including the aero constants) that produced the run, so a
snapshot is reproducible from source even though the data itself isn't committed.

## Regenerating

Run `testing/generateModelBaseline.m` (needs `physics/`, `visualization/`, and
`testing/helpers/` on the MATLAB path). It writes a fresh snapshot folder here.
