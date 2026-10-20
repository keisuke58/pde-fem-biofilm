# Real molar: the videos with the resolved stress (20 Oct 2026)

The three videos of the molar (`assets/anim_real_tooth_3d.mp4`, `anim_real_tooth_gate_3d.mp4`,
`anim_real_tooth_nutrient_gate_3d.mp4`) colour their third panel with the von Mises stress of the 32^3 full-box mesh
(h = 0.5 mm, the 0.5 mm film one element thick). `REAL_TOOTH_CONVERGENCE.md` shows that this stress is not converged (film
mean 2.0 mPa on that mesh against 3.3 mPa resolved; 95th percentile 5.3 against 10.5 mPa), and `REAL_TOOTH_SUBMODEL.md` /
`REAL_TOOTH_RESOLVED.md` give the converged values from the film submodel. `anim_real_tooth_resolved.py` renders the three
videos again with the von Mises panel from that submodel.

Script `ansys_usermat/apdl/anim_real_tooth_resolved.py`; videos `assets/anim_real_tooth_resolved_{base,gate,nutrient_gate}_3d.mp4`
with a last-frame `.png` each; numbers `real_tooth_videos.json` (this directory); the solved frames in
`figdata/real_tooth_resolved_<tag>_f<k>.npz` (ignored by git).

## What changed against the 32^3 videos

- The stress panel. Every third frame of each run (frames 0, 3, ..., 66; 23 per run, 69 in all) is solved again with
  `real_tooth_submodel.solve_submodel` on the run's own 96^3 field: rigid tooth, crop box around the film (the box of the
  final frame, 88 x 176 x 80 elements, so every frame of a run lives on one grid), void phi <= 0.05 dropped, wall = the
  field's voxel staircase, h = h_f/2 = 0.083 mm, six elements across the film, about 400k dof, 60 to 74 CG iterations,
  25 to 30 s per solve with the BLAS threads pinned to 1 and two solves at a time. Between the solved frames the stress
  is linear in time, as in the originals. Frame 0 has alpha - 1 = 0 everywhere and is stored as zero stress without a solve.
- Nothing else. The growth fields (phi, alpha - 1, c) are those of the run caches (`figdata/anim_real_tooth_3d.npz`,
  `real_tooth_gate.npz`, `real_tooth_nutrient_gate.npz`), so the first two panels are the same as before. The gate runs
  were driven by the 32^3 pressure during the run (p_h = 1.85 mPa, half the 32^3 film mean pressure of the base run);
  only their stress panel is re-solved, not their growth, as in `REAL_TOOTH_RESOLVED.md`. Camera, surface (phi > 0.05,
  smoothed), tooth, frames (67 + 24 turning), frame rate and the +30 degree end view of the PNG are those of the originals.
- The colouring of the film surface. As in the originals, each surface point shows the maximum of the stress along the
  inward normal over 1.5 mm (the stress through the film under that point; inactive cells, tooth and void, take the
  nearest film value, as the originals take the nearest film element). On the resolved mesh this maximum is the value of
  the first elements at the tooth wall, so the surface shows the near-wall stress (about 2.4 times the film mean,
  `REAL_TOOTH_SUBMODEL.md`) and the staircase of the voxel wall as fine stripes. The mean through the film instead
  (`--reduce mean`, weighted by the active elements on the line) was tried: median 2.2 mPa against 10.1 mPa for the
  maximum on the last frame of the base run, which on the common scale is a flat blue band with no structure, so the
  maximum is kept and stated on the frame note.
- The colour range: 0 to 14.2 mPa in all three videos, the 95th percentile of the value shown on the surface (the
  maximum through the film) on the last frame of the base run. The originals used the 95th percentile of the film
  elements of their own last frame (5.2, 2.2 and 0.6 mPa on the three 32^3 runs: each video had its own scale), and
  with a one-element film the element value and the shown value were the same thing. Here they are not: the
  volume-weighted 95th percentile over the film elements is 11.9 mPa (the converged number of `REAL_TOOTH_RESOLVED.md`),
  and with the maximum-through-the-film colouring a quarter of the surface lies above it, so the scale is taken on the
  shown value (5 % of the surface above the top). `--range elements` gives the 11.9 mPa scale; both numbers are in the
  JSON and on the frame note.

## Frames solved and wall time

| run | frames solved | solve time (23 solves) | render (91 frames) | film mean vM at the end (mPa) | band/core | 95th pct (mPa) | mean p (mPa) |
|---|---|---|---|---|---|---|---|
| base | 0, 3, ..., 66 (23) | 492 s | RENDER_BASE | 3.66 | 1.21 | 11.86 | 2.84 |
| gate | 0, 3, ..., 66 (23) | SOLVE_GATE | RENDER_GATE | GATE_LAST |
| nutrient + gate | 0, 3, ..., 66 (23) | SOLVE_NG | RENDER_NG | NG_LAST |

Two solve workers (`--solve --worker 0/1 --workers 2`, detached with nohup; each frame is written as soon as it is
solved and skipped on a rerun), WALL_SOLVE of wall time for the 69 frames in all; the renders one after the other,
WALL_RENDER. The end-of-run values reproduce `real_tooth_physics_resolved.json` to the last digit (same solve).

Over time (base run): film mean 0.64 mPa at T* = 0.05, 1.78 at 0.20, 2.59 at 0.40, 3.35 at 0.75, 3.66 at 1.10; the
95th percentile 1.3, 4.9, 8.2, 11.1, 11.9 mPa; band/core 1.11 rising to 1.29 at T* = 0.5 and back to 1.21. With the gate:
film mean 1.47 at 0.20, 1.98 at 0.40, 2.43 at 0.75, GATE_END at 1.10 (ratio to the base run RATIO_GATE at the end, as
`REAL_TOOTH_RESOLVED.md`). With nutrient and gate: NG_SERIES.

## Files

- `ansys_usermat/apdl/anim_real_tooth_resolved.py`: `--solve [--tags ...] [--worker k --workers n]` (frames striped over
  the workers), `--render [--tags ...] [--png-only] [--reduce max|mean] [--range surface|elements]`; without flags both,
  one solve at a time. Imports `real_tooth_submodel` (the solve) and `real_tooth_model` (geometry, caches, the surface
  colouring); the three original scripts are not changed.
- `assets/anim_real_tooth_resolved_{base,gate,nutrient_gate}_3d.mp4` and `.png` (last frame, +30 degree view as the
  original PNGs).
- `real_tooth_videos.json`: colour range and its two candidates, the solved frames and their times, the measures of every
  solved frame (film, band, core, layer means, percentiles, maximum, mean pressure) per run, render times.
- `figdata/real_tooth_resolved_<tag>_f<k>.npz`: vm, p (crop grid, nan outside the active elements), active, origin, h,
  shape, t, frame, measures and solver record (JSON strings), seconds.
