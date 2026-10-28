# Submodel of the film on the real molar: the stress across the film resolved (20 Oct 2026)

Question: `REAL_TOOTH_CONVERGENCE.md` found the growth field of the molar demonstration converged to 2 % but the
stress not converged on any full-box mesh (16^3 to 64^3, h = 1 to 0.25 mm): the film mean von Mises rose from 0.63 to
2.65 mPa with no end in sight, because the 0.5 mm film was half an element to two elements thick, and a full-box
mesh with five or more elements across the film (h <= 0.1 mm, 160^3 or finer) is beyond this machine. What are the
stress values of the model?

Script `ansys_usermat/apdl/real_tooth_submodel.py`; numbers in `real_tooth_submodel.json` (this directory); figure
`assets/fig_real_tooth_submodel.pdf` / `.png`; per-level fields (von Mises, p, masks, and the displacement for the
refined levels) in `figdata/real_tooth_submodel_<level>.npz`.

## What I did

Only the film is meshed, on the 192^3 field of the convergence study (h_f = 0.0833 mm; `real_tooth_conv_field_192.npz`),
with the mechanics of the demonstration (hex8 B-bar, small strain, E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain
(alpha - 1) I, pyamg smoothed aggregation + CG to a relative residual of 1e-10):

- Rigid tooth. The tooth (E = 1e4 Pa, a thousand times the film, clamped at the root) is replaced by u = 0 on every
  node of a tooth element. Checked on the full box at h = 0.333 mm against the elastic tooth: +0.10 % on the film
  mean, +0.26 % on the maximum.
- Crop. A box around the film: the bounding box of the field cells with phi > 0.01 at the final state plus a 1 mm
  margin, rounded to the field grid: x in [-8, -0.33], y in [-7.33, 7.67], z in [-5, 1] mm (7.7 x 15 x 6 mm). The
  background phi is k_alpha T* = 1.1e-3 everywhere, so 1e-3 is no threshold. Outside the film the material is the
  void of the demonstration, so the crop faces are free. Checked: crop against full box 0.00 %, 2 mm margin against
  1 mm 0.00 % (h = 0.333 and 0.167).
- Void dropped on the fine meshes. The void (E = 0.01 Pa against about 3 Pa in the film) is left out, the film
  surface is free: elements with phi <= 0.01 dropped change the film mean by 0.00 % (h = 0.333, 0.167, 0.083), phi <=
  0.05 by 0.02 % (h = 0.167; the mean pressure by 0.3 %). The levels with h >= h_f keep the void, the refined levels
  drop phi <= 0.05.
- Meshes aligned with the field grid: h = 4 h_f, 2 h_f (phi, alpha - 1 and the sets block-averaged), h_f (one
  element per cell) and h_f/2 (the field piecewise constant per cell). Trilinear interpolation from the cell
  centres instead of piecewise constant: -0.9 % on the film mean at h = h_f/2 on the 192^3 field (-3.4 % on the 96^3
  field at the same refinement ratio, consistent with a second-order field error).
- Wall. Two ways to classify the stress elements as tooth. Voxel wall: the field's own tooth (block majority above
  h_f, the containing cell below), i.e. the staircase the growth field was computed on, so from h = h_f on the
  geometry is fixed and only the mechanics is refined. STL wall: signed distance to the STL at the element centre,
  so the wall gets finer than the field's staircase below h_f (fine elements of a staircase tooth cell that lie
  outside the STL take the field of the nearest open cell). The two coincide at h = h_f. The STL wall moves the
  first-layer cells, which carry the highest stress, and changes the film mean by +2 % (h = 0.333), +4.6 %
  (0.167) and -6.2 % (0.042), the maximum by -21 % and the mean pressure by -12 % at 0.042. I take the voxel wall as
  the model (the field and the mechanics then see the same wall) and the STL-wall values as the sensitivity to the
  wall representation.
- A third point with fixed geometry, h = h_f/3 = 0.028 mm, on the whole film would need about 5.5M dof (the 2M-dof
  levels peak at 7.6 GB; pyamg's setup needs about 3.5 times the matrix), so it is solved on a sub-box around the
  middle of the flank (y in [-1.5, 1.5] mm, the whole z range of the film, void dropped) with the displacements of
  the h_f/2 solve prescribed on the two cut faces (trilinear from its nodes; the x and z faces are tooth or free),
  and the levels are compared on the same sub-region, the field cells with |y| < 0.7 mm (0.8 mm, 1.6 film
  thicknesses, from the cuts). Checks: the sub-box at the same h with the exact displacements reproduces the full
  solve to 1e-9; the sub-box at h_f/2 with the displacements of the h_f solve (a factor 2 coarser) differs from the
  full-film h_f/2 solve on the region by -0.13 % (film mean), -0.14 % (band mean), 0.00 % (maximum).
- Element sets as in `real_tooth_convergence.py`: film = field cells with phi > 0.3 that are not tooth, band = film
  with z in [-2.1, -1.1] mm (CEJ at -1.6), core = the rest; every stress element weighted by its fraction of field
  cells in the set, so the regions are the same on every mesh (film volume 17.33 mm^3 at h <= h_f). Added: the
  layer within 0.1 mm of the tooth surface (signed distance at the element centre), a physical region.
- A mesh-free reference from the field alone: for a flat film bonded to a rigid wall with a free surface the
  stress is local, von Mises = E (alpha - 1)/(1 - nu) and p = 2/3 of it in every cell. Its means: film 4.10 mPa, band
  4.43, core 3.80; p 2.74; 95th percentile 10.1 mPa.
- Code check: the script with the elastic tooth, voxel tooth and block means on the full 48^3 box reproduces the
  48^3 row of `real_tooth_convergence.json` to 1e-15.

## Convergence (whole film, voxel wall; von Mises in mPa, p = -tr(sigma)/3)

| level | h (mm) | elements across the film | dof | film mean | band mean | core mean | band/core | 0.1 mm at the tooth | 95th pct | 99th pct | maximum (at x, y, z mm) | mean p | max p | time |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 h_f | 0.333 | 1.5 | 36k | 2.443 | 2.564 | 2.328 | 1.10 | 1.84 | 6.58 | 8.51 | 9.79 (-4.2, 4.5, -2.2) | 2.15 | 29.8 | 4 s |
| 2 h_f | 0.167 | 3 | 277k | 2.762 | 2.916 | 2.617 | 1.11 | 4.51 | 8.59 | 10.1 | 12.6 (-3.9, 4.6, -2.4) | 2.04 | 42.9 | 28 s |
| h_f | 0.0833 | 6 | 2.17M | 3.257 | 3.464 | 3.062 | 1.13 | 7.33 | 10.40 | 12.6 | 15.7 (-4.8, 3.8, -2.2) | 2.75 | 119 | 276 s |
| h_f/2 | 0.0417 | 12 | 2.00M | 3.287 | 3.479 | 3.106 | 1.12 | 7.76 | 10.62 | 14.7 | 21.4 (-4.8, 3.8, -2.0) | 2.27 | 62 | 302 s |
| change per refinement | | | | +13.1 / +17.9 / +0.9 % | +13.7 / +18.8 / +0.4 % | +12.4 / +17.0 / +1.4 % | +1.1 / +1.5 / -1.0 % | +145 / +62 / +6 % | +31 / +21 / +2.1 % | +19 / +24 / +17 % | +28 / +25 / +36 % | -5 / +35 / -18 % | +44 / +178 / -48 % | |
| STL wall instead | | | | 2.495 / 2.889 / 3.257 / 3.081 | 2.615 / 3.041 / 3.464 / 3.266 | | | | | | 9.79 / 13.5 / 15.7 / 16.9 | 2.23 / 2.26 / 2.75 / 1.99 | | |
| full box, elastic tooth (`real_tooth_convergence.json`, 48^3 and 64^3) | 0.333, 0.25 | | | 2.440, 2.652 | 2.562, 2.752 | 2.326, 2.557 | | | 6.58, 7.43 | 8.51, 9.27 | 9.77, 11.4 | 2.15, 2.22 | 29.7, 38.3 | |

Middle of the flank (field cells with |y| < 0.7 mm; 1.854 mm^3 of film), fixed geometry, nested sub-boxes:

| level | h (mm) | dof | film mean | band mean | core mean | band/core | 0.1 mm at the tooth | 95th pct | 99th pct | maximum | mean p | max p | time |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h_f (from the whole-film solve) | 0.0833 | | 3.382 | 3.780 | 3.036 | 1.245 | 7.62 | 10.53 | 12.1 | 13.5 | 2.570 | 68 | |
| h_f/2 (from the whole-film solve) | 0.0417 | | 3.399 | 3.812 | 3.041 | 1.254 | 7.88 | 10.56 | 13.6 | 17.8 | 2.244 | 57 | |
| h_f/2, sub-box with the displacements of h_f (check) | 0.0417 | 425k | 3.395 | 3.806 | 3.038 | 1.253 | 7.88 | 10.56 | 13.6 | 17.8 | 2.242 | 57 | 43 s |
| h_f/3, sub-box with the displacements of h_f/2 | 0.0278 | 1.44M | 3.387 | 3.799 | 3.030 | 1.254 | 7.84 | 10.43 | 14.1 | 21.0 | 2.206 | 75 | 169 s |
| change per refinement | | | +0.5 / -0.35 % | +0.8 / -0.3 % | +0.2 / -0.4 % | +0.7 / 0.0 % | +3.4 / -0.5 % | +0.3 / -1.3 % | +12 / +4 % | +32 / +18 % | -12.7 / -1.7 % | -16 / +31 % | |
| flat-film estimate on the region | | | 4.17 | 4.73 | 3.69 | 1.28 | | | | | 2.78 | | |

What the tables say:

- The means converge once the film is resolved. From h_f on (six elements across the film, the field at its own
  resolution, the wall fixed) the film, band and core means change by +0.9, +0.4 and +1.4 % to h_f/2, and on the
  middle of the flank by +0.5 / -0.35 %, +0.8 / -0.3 % and +0.2 / -0.4 % over h_f, h_f/2, h_f/3. The changes are not
  monotone at the last step, so the three fine meshes resolve nothing beyond about 0.5 %; a Richardson fit is not
  meaningful and I quote the finest values with the spread of the last two steps as the error.
- The large steps before that are not the mechanics but the resolution of the field and the wall: from 2 h_f to
  h_f the film mean rises 18 %, because at 2 h_f the element holds the block mean of phi over eight cells (the mean of
  phi^2 is smaller) and the wall is the coarsened staircase. The free fit v0 + c h^p over the four whole-film levels
  (order 0.57, h -> 0 3.74 mPa) is therefore misleading, as was the "about 5 mPa" extrapolation of the full-box
  study: convergence is not a power law here; it is a change of regime at h = h_f. The full-box values (elastic
  tooth, block means) line up with the submodel at the same h: 2.440 (48^3) against 2.443 (4 h_f), 2.652 (64^3)
  between 4 h_f and 2 h_f.
- The 95th percentile (the colour scale of the video) is converged: 10.40, 10.62 mPa (whole film), 10.5, 10.6, 10.4
  (middle). The 99th percentile and the maximum are not: +17 % and +36 % at the last whole-film step, +12 / +4 % and
  +32 / +18 % on the middle, i.e. maximum ~ h^-0.4. The maximum sits on the tooth surface, at the +y end of the
  band region on the whole film ((-4.8, 3.8, -2.1) mm on every mesh, as in the full-box study) and in the band at
  the wall on the middle region; it is the stress concentration at the corners of the staircase wall, which is
  re-entrant at every step, so it grows without bound with h as the corner singularity of the cube study
  (`WALL_CONVERGENCE.md`). It is a value of the staircase, not of the model.
- The mean pressure converges from above: 2.75 (h_f), 2.27 (h_f/2) on the whole film; 2.57, 2.24, 2.21 on the
  middle (order 2.8, h -> 0 2.19 mPa). At h_f the first element in each notch of the staircase has most of its nodes
  fixed and carries an almost hydrostatic stress, which inflates the mean by 20 %; one refinement removes most of
  it. The maximum pressure is a corner value like the maximum von Mises.
- The 0.1 mm next to the tooth: 7.33, 7.76 mPa (whole film), 7.62, 7.88, 7.84 (middle), 2.3 to 2.4 times the film
  mean; converged to about 1 % from h_f/2 on (at h = 0.333 and 0.167 the layer is one element, which is why it rises
  by 145 and 62 % there).
- Band/core: 1.13, 1.12 on the whole film, 1.25 on the middle of the flank (where the band is the middle of the
  seed in z and the ends of the seed are outside the region): the band is 12 % above the rest of the film overall,
  25 % in the middle of the flank.
- Against the flat-film estimate the FEM means are 20 % lower (film 3.29 against 4.10; middle 3.39 against 4.17; band
  3.48 against 4.43), the 95th percentile 5 % higher and the mean pressure 17 % lower. The difference is the finite
  film: it is 2.4 mm long and 0.5 mm thick, so the end zones where the biaxial stress relaxes are a large part of it,
  the flank is curved (t/R about 0.1) and the wall is a staircase; the two wall representations alone differ by 6 %.
  The flat-film number is an upper bound that checks the order of magnitude, not a reference value.
- The field grid at fixed h = 0.083 mm (96^3 against 192^3 field, both with void dropped): film mean -0.9 %, band
  +2.4 %, core -5.3 %, mean pressure -16 %, maximum +0.6 %. So for the comparisons between runs of
  `real_tooth_model.py` (96^3 field) the film mean at h = 0.083 mm is good to 1 %, the band and core means to about
  5 %, and the mean pressure should be compared only between runs on the same grid.

## Converged values

Model: the film of the demonstration on the staircase wall of its field (192^3 grid, 0.083 mm), rigid tooth. The error
is the spread over the last two refinements (mesh) and, in brackets, the change when the wall is taken from the STL
instead of the field's voxels (geometry representation):

- film mean von Mises 3.3 mPa (+- 1 %; wall -6 %), band mean 3.5 mPa (+- 1 %; -6 %), core mean 3.1 mPa (+- 1.5 %;
  -6 %), band/core 1.12 (+- 1 %); middle of the flank: film 3.39, band 3.80, core 3.03 mPa (+- 0.5 %), band/core 1.25;
- 95th percentile 10.5 mPa (+- 2 %; wall -10 %); the 0.1 mm next to the tooth 7.8 mPa (+- 1 %; -7 %);
- mean pressure 2.2 to 2.3 mPa (the h_f/2 value 2.27 with about 3 % left, the middle-region limit 2.19; wall -12 %);
- the maximum von Mises and the maximum pressure have no converged value (~h^-0.4; 21 mPa and 62 to 75 mPa at the
  finest meshes); the 99th percentile neither (14 to 15 mPa, +4 % at the last step).

What this means for the demonstration (`REAL_TOOTH.md`, the video on the 32^3 mesh): its film mean of 2.0 mPa is
40 % below the value of the model, 3.3 mPa; its maximum of 6 mPa is the one-element value of a quantity that has no
limit; its 95th percentile of 5.3 mPa is half the converged 10.5 mPa. The statements that survive: the film stress
is of the order of a few mPa (3.3 mPa film mean), roughly uniform along the film with the cervical band 12 % above
the rest (25 % in the middle of the flank), 2.4 times higher in the 0.1 mm next to the tooth than on average, and
with its largest values at the tooth surface at the ends of the band region on the flank. The comparison with the
2 mm cube (stresses of the same order, mPa) holds.

## Timing and memory (this machine, four cores, BLAS threads 1, the field solves of other jobs running alongside)

| level | dof | prepare | assembly | AMG setup | CG | total | peak RSS |
|---|---|---|---|---|---|---|---|
| 4 h_f (0.333 mm) | 36k | 0 s | 1 s | 1 s | 2 s | 4 s | |
| 2 h_f (0.167 mm) | 277k | 0 s | 6 s | 8 s | 14 s | 28 s | |
| h_f (0.083 mm, void kept) | 2.17M | 0 s | 56 s | 90 s | 130 s | 277 s | 7.6 GB |
| h_f/2 (0.042 mm, void dropped) | 2.00M | 1 s | 61 s | 90 s | 147 s | 302 s | about 7 GB |
| h_f/3 on the sub-box | 1.44M | 0 s | 35 s | 57 s | 75 s | 169 s | about 5 GB |
| h_f/2 on the 96^3 field (0.083 mm) | 402k | 14 s | 7 s | 10 s | 24 s | 55 s | |

The STL wall costs the signed-distance evaluation at the element centres (93 s at 9.5M elements). pyamg's default
dense coarse solver (pinv) took 20 to 90 s on the small systems whose coarsest level has a few thousand unknowns;
the script uses a sparse LU there (same result, 3 s). 46 to 68 CG iterations at every size.

## Files

- `ansys_usermat/apdl/real_tooth_submodel.py`: `load_field(192 | "gate" | path)` reads a convergence field cache or
  a run cache of `real_tooth_model.py` (final frame, 96^3) or `anim_real_tooth_3d.py` ("base"); `solve_submodel(phi,
  am1, h, tooth=, sdf=, h_f=, box=, margin=, void=, void_thr=, interp=, tooth_mode=, rigid=, u_bc=, region=)` returns
  the measures, the solver record with the times, and the fields on the crop grid. `--checks`, `--levels m4 m2 r1
  r2`, `--subbox r2b r3`, `--refit`, `--figure`; every solve goes to the JSON and an npz as soon as it is finished and
  is skipped on a rerun. For the gate and nutrient runs (96^3 field): `solve_submodel(..., h=L/96/2)` (h = 0.083 mm,
  about 400k dof, one minute) gives the film mean to 1 % of the 192^3-field value.
- `real_tooth_submodel.json`: model, field (crop box, flat-film reference), checks with labels, the checks summary
  (percent differences), the levels (whole film, both walls), the fits per sequence, the sub-box levels with the
  cut-boundary check.
- `assets/fig_real_tooth_submodel.pdf`, `.png`: (a) the measures against h (log-log): submodel with the voxel wall
  solid, STL wall dashed, the full-box values of `real_tooth_convergence.py` dotted, the middle-of-the-flank sequence
  dash-dotted with stars, the flat-film film mean as a horizontal line; (b) the von Mises stress on the h_f/2 mesh in
  the section z = -1.6 mm through the band, tooth grey; (c) the vertical section y = 0 across the film thickness.
- `figdata/real_tooth_submodel_<level>.npz` (ignored by git): vm, p, tooth, active, w_film, sdf, origin, h, shape, and
  u for r1, r2, r3 (the sources of the cut boundaries).
