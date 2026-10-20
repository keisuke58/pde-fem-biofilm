# Stress where the seed meets the clamped face: what converges with the mesh (20 Oct 2026)

Question: the 3D demonstrations with the seed in a corner or on an edge of the cube and the bottom face clamped
(`anim_seed_offcentre_3d.py`, slides `slides_3d_demo.tex`) show their largest von Mises stress at the clamped edge
under the seed, about 15 mPa on the 16^3 stress mesh against 0.76 mPa for the seed in the centre. Is that a value of
the model or a value of the mesh?

Script `ansys_usermat/apdl/wall_stress_convergence.py`; numbers in `wall_convergence.json` (this directory);
figure `assets/fig_wall_convergence.pdf` / `.png`.

## What I did

- One growth field per case on a 96^3 grid (RK4, the solver of `wp3_reference.py`: beta = 0.02 mm^2/T*,
  k_alpha = 1e-3 /T*, penalty 5, zero-flux walls, T* = 1.1). 96 is a multiple of 16, 24, 32 and 48, so phi and
  alpha - 1 are block means of the same field on every stress mesh; only the mechanics changes between the meshes.
  (The demos used a 48^3 field, which 32^3 does not divide; the 16^3 stresses here agree with the demos to 0.1 %.)
- Mechanics on 16^3, 24^3, 32^3 and 48^3: hex8 with B-bar, E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain
  (alpha - 1) I, small strain. Corner and edge seed: bottom face z = -1 clamped. Centre seed (the Chapter 4 problem,
  reference): only the three corner-node rigid-body constraints of the deck. Direct solve up to 24^3, Jacobi-
  preconditioned CG on 32^3 and 48^3 (relative residual 1e-10, 1200 to 1900 iterations with the clamped face, 3800
  and 5900 with the rigid-body constraints; 48^3 takes 11 to 15 min).
- Measures (biofilm = elements with phi > 0.5; p = -tr(sigma)/3, compression positive):
  seed mean (initial seed region); band = biofilm with z in [-1, -0.75], the quarter millimetre on the clamped face;
  core = biofilm with z in [-0.75, -0.5]; wall band = biofilm within 0.25 mm of a free wall the seed touches, above
  the band; maximum with its location; 95th and 99th percentiles over the biofilm; maximum p. Same bands on every
  mesh, as in the cervical study of `p3/uel/README.md`.

## Results (von Mises in mPa; the mesh column is the stress mesh)

Corner seed (touches the clamped face and the walls x = -1, y = -1):

| mesh | seed mean | band mean | core mean | band/core | wall-band max | 95th pct | 99th pct | maximum | max p |
|---|---|---|---|---|---|---|---|---|---|
| 16^3 | 4.42 | 7.76 | 1.45 | 5.35 | 3.37 | 16.0 | 16.2 | 16.5 | 9.6 |
| 24^3 | 4.37 | 7.64 | 1.46 | 5.23 | 3.28 | 17.4 | 18.9 | 19.5 | 12.7 |
| 32^3 | 4.35 | 7.60 | 1.48 | 5.14 | 4.29 | 17.6 | 20.6 | 21.8 | 15.4 |
| 48^3 | 4.33 | 7.61 | 1.49 | 5.11 | 4.83 | 16.7 | 21.8 | 26.3 | 20.1 |
| h -> 0 | 4.29 (order 1.1) | 7.5 to 7.6 | 1.51 (p = 1) | 5.03 (order 1.4) | rises, ~h^-0.4 | no trend | ~h^-0.27 | ~h^-0.42 | ~h^-0.67 |

Edge seed (touches the clamped face and the wall x = -1):

| mesh | seed mean | band mean | core mean | band/core | wall-band max | 95th pct | 99th pct | maximum | max p |
|---|---|---|---|---|---|---|---|---|---|
| 16^3 | 3.59 | 6.40 | 1.31 | 4.90 | 1.99 | 13.7 | 15.0 | 15.1 | 8.8 |
| 24^3 | 3.57 | 6.32 | 1.32 | 4.79 | 2.11 | 15.0 | 17.0 | 17.9 | 10.3 |
| 32^3 | 3.56 | 6.34 | 1.33 | 4.76 | 2.70 | 14.9 | 17.8 | 19.9 | 12.2 |
| 48^3 | 3.55 | 6.36 | 1.34 | 4.75 | 3.13 | 14.3 | 18.5 | 22.4 | 15.4 |
| h -> 0 | 3.54 (order 1.1) | 6.3 to 6.4 | 1.39 (order 0.5) | 4.74 (order 2.9) | rises, ~h^-0.4 | no trend | ~h^-0.19 | ~h^-0.36 | ~h^-0.52 |

Centre seed, rigid-body constraints only (the reference):

| mesh | seed mean | 95th pct | maximum | max p |
|---|---|---|---|---|
| 16^3 | 0.618 | 0.76 | 0.76 | 1.42 |
| 24^3 | 0.615 | 0.74 | 0.79 | 1.51 |
| 32^3 | 0.613 | 0.75 | 0.79 | 1.54 |
| 48^3 | 0.612 | 0.75 | 0.79 | 1.57 |
| h -> 0 | 0.609 (order 1.2) | converged | converged (-0.2 % at the last step) | +1.4 % at the last step |

The seed mean of the centre case agrees with `wp3_stress_reference.json` (0.6176, 0.6143, 0.6113 mPa on 16^3,
24^3, 48^3 there) to 0.1 %, so the 96^3 field is as good as the Richardson reference field for this purpose.

## What converges and what does not

Converging (changes between 32^3 and 48^3 below 1 %):

- The seed mean (order about 1.1, as in Chapter 4): corner 4.3 mPa, edge 3.5 mPa, centre 0.61 mPa.
- The band mean, the mean over the bottom quarter millimetre of the biofilm: corner 7.6 mPa, edge 6.4 mPa
  (16^3 is 2 % and 1 % above the 48^3 value; from 24^3 on the changes are within 0.6 % and not monotone, so the
  four meshes resolve nothing beyond "7.6" and "6.4"). The band mean of p behaves the same (4.1 and 3.7 mPa).
- The core mean (the quarter millimetre above the band): 1.5 and 1.3 mPa, rising by 0.6 to 1 % per refinement,
  order about 0.5; the first-order extrapolation is 2 to 4 % above the 48^3 value.
- The ratio band mean / core mean: 5.1 (corner) and 4.75 (edge), extrapolated 5.0 and 4.7. "About 5" is safe.
- The 95th percentile of the von Mises stress over the biofilm, which is what the colour scale of the videos
  shows: corner 16.0, 17.4, 17.6, 16.7 mPa; edge 13.7, 15.0, 14.9, 14.3 mPa. No trend after 24^3, scatter of
  about 5 %. With the demo's own threshold phi > 0.3 the same percentile is 15.4, 16.3, 15.5, 15.1 (corner) and
  13.0, 13.1, 12.2, 12.5 mPa (edge). The 95th percentile of p (the base of p_h in the growth-law demo) is also
  stable: 8.5 to 9.2 mPa (corner), 7.7 to 8.2 mPa (edge).
- The centre case altogether: its maximum, 0.76 mPa on 16^3, is 0.79 mPa from 24^3 on.

Not converging:

- The maximum von Mises stress: corner 16.5, 19.5, 21.8, 26.3 mPa, edge 15.1, 17.9, 19.9, 22.4 mPa, growing by
  12 to 20 % per refinement, ~h^-0.42 (corner) and ~h^-0.36 (edge) over the four meshes, with no sign of
  levelling off (the pair slopes are 0.41, 0.40, 0.45 for the corner). The maximum of p grows faster, ~h^-0.67 and
  ~h^-0.52. The 99th percentile follows the maximum at a smaller rate (+6 % and +4 % at the last refinement).
- Where it sits: not at the boundary of the seed but at the corner of the box. Corner case: in the element at the
  vertex (-1, -1, -1) from 32^3 on (the first elements at 16^3 and 24^3 are the same place at their resolution),
  where the clamped face meets two free walls. Edge case: on the line x = -1, z = -1, where the clamped face meets
  the free wall, at y about 0, the middle of the seed. This is a corner of the elastic problem (a clamped face
  meeting a free face at a right angle inside a material that is pushed outward by growth), and the element value
  next to it grows like a power of the distance h/2 to the corner. The exponent of the corresponding wedge problem
  depends on nu and I have not computed it; the observed 0.4 is what one expects for such a corner, not an error
  of the solver (CG residual 1e-10; the direct solves on 16^3 and 24^3 show the same growth).
- The maximum of the wall band (the free re-entrant edges above the band): 3.4 to 4.8 mPa (corner), 2.0 to
  3.1 mPa (edge), also at the wall and in the first element above z = -0.75, i.e. the tail of the same
  clamped-edge rise running up the wall. It grows at about the same rate but is five times smaller than the
  maximum at the clamped edge; the wall band's mean (1.7 and 1.6 mPa) converges. The free walls alone do not
  produce a concentration worth a number on the slides.

So the picture of the cervical study holds here: a fixed physical band next to the constrained face has a
converged mean, the pointwise maximum at the corner of the elastic problem does not, and a ratio stated with the
maximum is a ratio stated with the mesh.

## Ratios that can be stated

With the converged measures (48^3 values; the extrapolations change them by at most 2 %):

- band mean / centre seed mean: corner 7.61 / 0.612 = 12.4, edge 6.36 / 0.612 = 10.4: "10 to 12 times".
- band mean / core mean: 5.1 and 4.75: "about 5 times the layer above".
- seed mean / centre seed mean: 7.1 (corner), 5.8 (edge).
- 95th percentile: about 17 mPa (corner), 14 mPa (edge), "about 15 mPa" on the colour scale, mesh-independent
  within 5 %.

The ratios in the slides come from the 16^3 maxima: 16.5 / 0.76 = 22 ("about 20 times the centre value") and
"about 15 mPa" as the largest value. On 48^3 the same ratio is 26.3 / 0.79 = 33 and keeps growing, so neither
number is a property of the model.

## What to change in `slides_3d_demo.tex`

Frame "2. Seed in a corner and on an edge", now: "the stress concentrates where the seed meets the clamped face,
about 15 mPa against 0.76 mPa in the centre. The largest values sit at the clamped edge and depend on the mesh; the
colour scale stops at the 95th percentile."

Suggested: "The stress concentrates in the biofilm next to the clamped face: the mean over its bottom quarter
millimetre is 7.6 mPa (corner) and 6.4 mPa (edge), 10 to 12 times the seed stress of the centre case (0.61 mPa)
and about 5 times the mean of the layer above it. These means change by less than 1 % between the 24^3 and 48^3
stress meshes. The largest values sit where the clamped face meets a free wall and grow with refinement (16 to
26 mPa from 16^3 to 48^3, about h^-0.4), as at a corner of an elastic problem; the colour scale stops at the 95th
percentile, about 15 mPa, which does not depend on the mesh."

Frame "3. ... growth law", now: "Without the law the stress at the clamped edge reaches about 15 mPa; with the law
it stays well below." Suggested: "Without the law the stress in the biofilm next to the clamped face reaches about
15 mPa (95th percentile, the top of the colour scale); with the law it stays well below." (p_h = 3.3 mPa was set
from the 95th percentile of p on 16^3, 8.5 mPa; on 48^3 that percentile is 8.7 mPa, so p_h is not a mesh value.)

Summary frame, now: "the stress then concentrates where growing biofilm meets a fixed surface, about 20 times the
centre value." Suggested: "the stress then concentrates in the biofilm next to the fixed surface, 10 to 12 times
the seed stress of the centre case (band mean over the bottom quarter millimetre; the pointwise maximum at the
corner of the box is mesh-dependent and is not quoted)."

Numbers that can be kept as they are: "about 15 mPa" when it is said to be the 95th percentile or the colour scale,
0.76 mPa only as "the largest value of the centre case on 16^3" (0.79 mPa from 24^3 on; the seed mean 0.61 mPa is
the better reference), and the 16^3 / 48^3 mesh labels in the table. Numbers to drop: "about 20 times" and any
maximum at the clamped edge given as a value of the model.

## Files

- `ansys_usermat/apdl/wall_stress_convergence.py` (solve per case, `--refit`, `--figure`).
- `ansys_usermat/apdl/results/2026-10-3d/wall_convergence.json`: every measure per case and mesh, the solver
  record (iterations, residual), the fits (change per refinement, first- and second-order extrapolations, free
  fit with its order and a reliability flag, power-law exponent and pair slopes for the maxima).
- `assets/fig_wall_convergence.pdf`, `.png`: (a) corner and (b) edge, the measures against h on log axes; (c) the
  ratios to the centre seed mean and band/core.
- Caches in `figdata/`: `wall_field_<case>_96.npz` (the fields), `wall_stress_<case>_<n>.npz` (von Mises, p and
  phi per element), for further post-processing without a new solve.
