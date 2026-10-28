# Mesh convergence of the real-molar demonstration (20 Oct 2026)

Question: the demonstration on the molar FDI 46 (`REAL_TOOTH.md`, `anim_real_tooth_3d.py`) computes its growth
field on 96^3 cells (h = 0.167 mm) and its stress on 32^3 hex8 elements (h = 0.5 mm, the film one element
thick). Which of its numbers are values of the model and which are values of the grid?

Script `ansys_usermat/apdl/real_tooth_convergence.py`; numbers in `real_tooth_convergence.json` (this
directory); figure `assets/fig_real_tooth_convergence.pdf` / `.png`; caches `figdata/real_tooth_conv_field_<n>.npz`
(the fields) and `figdata/real_tooth_conv_stress_f<field>_m<mesh>.npz` (von Mises, p, phi per element).

Status at the time of writing: the field stage is complete (48^3, 96^3, 144^3, 192^3, plus 96^3 with dt/2), the
stress stage is complete on the 96^3 field (16^3, 24^3, 32^3, 48^3) and on the 192^3 field up to 48^3; the
64^3 mesh on the 192^3 field (h = 0.25 mm, 823k dof) was still running (it is written to the JSON when it finishes). To finish after the chain:
`python ansys_usermat/apdl/real_tooth_convergence.py --fields 48 96 144 192 --meshes 16 24 32 48 64
--stress-fields 96 192` (everything comes from the caches and the JSON, only missing pairs are solved; it also
writes the dt/2 check and the field fits into the JSON), then `--figure`.

## What I did

- Field: the model and geometry of the demonstration (beta = 0.08 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5,
  T* = 1.1, RK4, zero flux at the walls and at the voxelised tooth, seed phi = 1 within 0.5 mm of the tooth on the
  -x flank at z = -3 .. -1 mm), final state on 48^3, 96^3, 144^3 and 192^3 cells (h = 0.333, 0.167, 0.111,
  0.083 mm). The tooth and the seed are re-voxelised on every grid, so the field convergence includes the
  staircase geometry. dt = min(0.15 h^2/beta, T/66); 96^3 was repeated with dt/2.
- Stress: the mechanics of the demonstration (hex8 B-bar, E = 10 Pa (phi^2 + 1e-3), tooth E = 1e4 Pa without
  growth, nu = 0.49, eigenstrain (alpha - 1) I, bottom face clamped; pyamg smoothed aggregation + CG to a
  relative residual of 1e-10, 62 to 82 iterations) on 16^3, 24^3, 32^3 and 48^3 elements (h = 1, 0.667, 0.5,
  0.333 mm) with phi and alpha - 1 block-averaged from the 96^3 field, and on 16^3 to 64^3 from the 192^3 field.
  Every stress mesh divides the field grid, so between two meshes only the mechanics changes. 32^3 on the 96^3
  field is the demonstration itself (same maximum, 6.0 mPa, in the same element).
- Element sets in physical space, the same on every mesh. They are defined on the field grid (film = cells with
  phi > 0.3 that are not tooth; cervical band = film with z in [-2.1, -1.1] mm, 0.5 mm either side of the CEJ at
  z = -1.6; core = the rest of the film) and each stress element is weighted by the fraction of its field cells in
  the set. The means are then volume means over the same region on every mesh. With the demonstration's own
  definition (block-averaged phi > 0.3) the 16^3 mesh would hold 13 film elements and no band, because a 1 mm
  element with half a millimetre of film in it averages below 0.3. Maximum, its location and the element counts
  use the majority-film elements (fraction > 0.5); the percentiles are volume-weighted over every element that
  holds film. First layer = majority-film elements next to a tooth element, outer = the rest; that layer is one
  element thick, so it is a different region on every mesh and I report it only for completeness. Elements that
  are majority tooth are tooth, so the film inside them is lost on the coarse meshes (film volume 13.5 mm^3 on
  16^3 against 16.6 mm^3 on 48^3 and 17.2 mm^3 on the field grid).

## Field convergence

| grid | h (mm) | film volume (sum phi h^3, mm^3) | volume phi > 0.5 (mm^3) | seed mean phi | seed mean alpha - 1 | max alpha - 1 | seed volume (mm^3) |
|---|---|---|---|---|---|---|---|
| 48^3 | 0.333 | 17.120 | 9.26 | 0.566 | 7.970e-4 | 1.0324e-3 | 13.56 |
| 96^3 | 0.167 | 17.241 | 8.61 | 0.542 | 7.569e-4 | 1.0142e-3 | 13.68 |
| 144^3 | 0.111 | 17.238 | 8.33 | 0.538 | 7.499e-4 | 1.0068e-3 | 13.67 |
| 192^3 | 0.083 | 17.22 | 8.25 | 0.536 | 7.468e-4 | 1.005e-3 | 13.66 |
| 96^3, dt/2 | 0.167 | 17.24 | 8.61 | 0.542 | 7.569e-4 | 1.014e-3 | 13.68 |

- The film volume (the integral of phi, the quantity the growth law conserves) is converged: 96^3, 144^3 and
  192^3 agree within 0.1 %; 48^3 is 0.7 % low. The seed and tooth volumes (staircase geometry) are within 0.1 %
  from 96^3 on.
- The seed mean of alpha - 1 changes by -5.0 %, -0.9 % and -0.4 % per refinement; the Richardson order on
  48/96/192 is 2.0 and the h -> 0 value 7.43e-4, so the 96^3 value of the demonstration (7.57e-4) is 1.8 % high.
  The seed mean phi behaves the same (order 2.1, h -> 0 0.535, 96^3 1.3 % high).
- The maximum of alpha - 1 converges at first order (1.032, 1.014, 1.005e-3; Richardson order 1.0, h -> 0
  0.996e-3; 96^3 1.9 % high). It sits at the tooth surface, where the staircase changes with the grid.
- The volume with phi > 0.5 (a threshold on the diffuse front) still changes by 3 % and 1 % at the last two
  refinements (order about 0.5); it is a value of the grid at this resolution and I do not quote it.
- The time step is not an issue: dt/2 on 96^3 changes every measure by less than 0.01 %.

So the growth field of the demonstration is converged to about 2 % in its means and in its maximum, and the
96^3 grid is adequate for the field. The stress, however, is not limited by the field but by the stress mesh.

## Stress convergence (96^3 field; von Mises in mPa, p = -tr(sigma)/3)

| mesh | h (mm) | film elements | film mean | band mean | core mean | band/core | first layer | outer | 95th pct | 99th pct | maximum (at x, y, z mm) | max p | mean p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16^3 | 1.0 | 11 | 0.66 | 0.60 | 0.72 | 0.84 | 1.04 | - | 1.73 | 1.77 | 1.77 (-4.5, 4.5, -2.5) | 3.7 | 0.85 |
| 24^3 | 0.667 | 43 | 1.28 | 1.42 | 1.14 | 1.25 | 1.77 | 1.14 | 3.63 | 4.22 | 4.23 (-4.33, 4.33, -1.67) | 11.3 | 1.38 |
| 32^3 | 0.5 | 133 | 2.05 | 2.17 | 1.93 | 1.13 | 2.51 | 1.30 | 5.26 | 5.96 | 6.00 (-3.75, 4.75, -2.25) | 16.5 | 2.03 |
| 48^3 | 0.333 | 389 | 2.43 | 2.62 | 2.24 | 1.17 | 3.97 | 2.02 | 6.64 | 8.42 | 10.14 (-4.17, 4.5, -2.17) | 34.1 | 2.03 |
| change per refinement | | | +93 / +61 / +18 % | +136 / +53 / +21 % | +57 / +70 / +16 % | +50 / -10 / +4 % | +69 / +42 / +58 % | | +109 / +45 / +26 % | +139 / +41 / +41 % | +139 / +42 / +69 % | +208 / +46 / +106 % | +63 / +47 / +0.1 % |
| fit over the four meshes | | | order 0.45, h -> 0 5.4 | order 0.66, h -> 0 4.6 | not monotone in order | | | | ~h^-1.2 | ~h^-1.4 | ~h^-1.6 | ~h^-2.0 | |

On the 192^3 field the same meshes give film means 0.63, 1.23, 1.93, 2.44 mPa, band means 0.57, 1.31, 2.01,
2.56 mPa and maxima 1.73, 4.09, 5.78, 9.77 mPa: within 5 % of the 96^3-field values on every mesh and within
0.5 % for the film mean on 48^3, so the field grid is not what limits the stress; the stress mesh is. The 64^3
mesh on the 192^3 field (h = 0.25 mm, two elements across the film, 811k dof, 117 CG iterations, 100 s) gives
film mean 2.65 mPa (+8.7 % on 48^3), band 2.75 (+7.4 %), core 2.56, band/core 1.08, 95th percentile 7.43,
99th 9.27, maximum 11.4 mPa at (-4.9, 3.9, -2.1) mm (+16.5 %, pair slope 0.53 against 1.3 before), maximum p
38 mPa, mean p 2.22 mPa. Over the five meshes the film mean has observed order 0.51 and a free-fit h -> 0 value
of 4.8 mPa, the band mean order 0.66 and 4.3 mPa; the last changes (7 to 9 %) are the first below 10 %. The
growth of the maximum slows once the film is two elements thick (~h^-1.35 over the five meshes). The next step
is the 96^3 mesh (h = 0.167 mm, 2.7M dof), which the chunked assembly makes possible on this machine.

What the table says:

- None of the stress measures is converged on the meshes of the demonstration. The film mean rises by 93, 61
  and 18 % per refinement, the band and core means likewise; the free fit gives an observed order of 0.45 to
  0.66 and an h -> 0 value about twice the 48^3 value, but with changes of this size the fit is an
  extrapolation of a trend, not a converged number. The reason is the film thickness: 0.5 mm, which is half an
  element on 16^3, one element on 32^3 and one and a half on 48^3. A thin layer under an eigenstrain carries its
  stress in the gradient across its thickness, and an element that holds the whole layer (plus space or tooth)
  averages that gradient away. The cube study (`WALL_CONVERGENCE.md`), whose film is 0.25 mm on a 2 mm cube,
  has 2 to 6 elements across the film on the same mesh counts, and its band means converge to 1 %; the molar
  at real size needs the same number of elements across 0.5 mm, i.e. h of 0.1 mm or less (160^3 on this box),
  which pyamg can do on this machine only up to about 96^3 (2.7M dof).
- The maximum grows faster than in the cube, ~h^-1.6 against h^-0.4 there, and the maximum of p ~h^-2. This is
  not a corner singularity of the elastic problem: the growth of the maximum follows the under-resolved film
  (the first-layer mean grows at the same rate, +58 % at the last step), and I expect it to slow once the film is
  resolved. The location of the maximum is stable in physical space: the end of the band on the +y side at
  (-4.2 +- 0.4, 4.5 +- 0.3, -2.2 +- 0.3) mm on every mesh, where the film ends on the curved flank. The 95th
  percentile (the colour scale of the video) rises by 26 % at the last refinement.
- The ratio band / core is 1.1 to 1.25 from 24^3 on (not monotone, +4 % at the last step): the cervical band
  is not a concentration in this demonstration (the seed covers z = -3 .. -1, the band z = -2.1 .. -1.1, so band
  and core are both inside the seed, unlike the cube where the band sits on the clamped face). "The stress is
  roughly uniform along the film, 10 to 25 % higher in the band" is as much as the ratio supports.
- The mean pressure in the film is the one measure that stopped moving: 2.03 mPa on 32^3 and 48^3 (+0.1 %);
  its band mean changed by 4 %.

## What is safe to state about the real-tooth result

- The growth field: film volume 17.2 mm^3 (converged to 0.1 %), mean alpha - 1 in the seed 7.5e-4 (96^3 value
  1.8 % high, Richardson order 2), maximum alpha - 1 1.0e-3 (2 % high, order 1). These are values of the
  model.
- The stress level: "of the order of millipascal", "the film mean is a few mPa", "the stress is roughly uniform
  along the film with its largest values at the ends of the band". The film mean von Mises on the meshes
  computed here is 0.66, 1.28, 2.05 and 2.43 mPa (16^3 to 48^3) and extrapolates to about 5 mPa, so a sentence
  needs the mesh with it: "2.0 mPa on the 32^3 mesh of the video, 2.4 mPa on 48^3, 2.65 mPa on 64^3, not
  converged (observed order about 0.5, still +9 % at the last step; extrapolated 4 to 5 mPa)". The maximum of
  6.0 mPa in `REAL_TOOTH.md` and in the video is the 32^3 value; on 48^3 it is 10.1 mPa and on 64^3 11.4 mPa,
  and it should not be quoted as a value of the model. The same holds for the maximum pressure (16.5 mPa on
  32^3, 34 to 38 mPa on 48^3 and 64^3) and for the 95th percentile (5.3, 6.6, 7.4 mPa).
- The comparison with the 2 mm cube ("stresses of the same order, mPa") holds; a number beyond the order does
  not.
- The location of the maximum, the end of the band on the +y side, is mesh-independent.

I would therefore word the demonstration as: "The von Mises stress in the film is of the order of a few mPa,
roughly uniform along the band (the cervical band is 10 to 25 % above the rest), with its largest values at the
ends of the film on the curved flank. The stress mesh of the video (0.5 mm, one element across the film) does not
resolve the stress across the film; the film mean rises from 2.0 mPa on that mesh to 2.4 mPa on 48^3 and is not
converged, so the maximum (6 mPa on the video's mesh) is a value of the mesh, not of the model. The growth field
itself is converged to 2 %." And in `REAL_TOOTH.md`, the limitation "the stress mesh is coarse" should say that
the stress numbers are not converged, with the numbers above.

## Files

- `ansys_usermat/apdl/real_tooth_convergence.py`: `--field-only`, `--stress-only` (fields from the caches),
  `--meshes`, `--stress-fields`, `--refit`, `--figure`, `--force-stress`. Each (field, mesh) pair is written to the
  JSON and to its npz as soon as it is solved and skipped on a rerun, so a restart loses at most one solve. The
  stiffness matrix is assembled in element chunks straight into the free-dof matrix, so a 96^3 stress mesh
  (2.7M dof) fits in memory.
- `results/2026-10-3d/real_tooth_convergence.json`: model, field measures per grid, every stress measure per
  field and mesh with the solver record, the fits (change per refinement, first- and second-order extrapolations,
  free fit with order and reliability flag, power-law exponent and pair slopes for the maxima and percentiles,
  field influence at fixed mesh).
- `assets/fig_real_tooth_convergence.pdf`, `.png`: (a) the measures against the element size on log axes,
  96^3 field solid, 192^3 field dashed; (b) the ratios (maximum, 99th and 95th percentile to the film mean,
  band/core, first layer/outer).
