# Growth stress in the measured 3D structure of the biofilm (CLSM stacks, 20 Oct 2026)

Question: what does the real, heterogeneous structure of the biofilm do to the growth stress, against a
laterally uniform film of the same thickness profile? Ten CLSM/FISH stacks (commensal CH and dysbiotic DH,
days 1, 6, 10, 15, 21) are placed on a flat, clamped substratum, grown with the thesis growth model and
solved for the stress on the voxel mesh.

Script `ansys_usermat/apdl/clsm_microscale.py` (`run`, `checks`, `figure`, `video`, `table`); numbers in
`clsm_microscale.json` (this directory); figure `assets/fig_clsm_microscale.pdf` / `.png`; video
`assets/anim_clsm_microscale_DH_d21_3d.mp4` / `.png`; computed fields (von Mises, pressure, alpha - 1) in
`figdata/clsm_microscale_*.npz` (not tracked).

## Data source and processing

The stacks are the five-species in-vitro oral biofilms of N. Heine et al. (MHH/NIFE, HOBIC flow chamber,
titanium substrate), imaged by CLSM after FISH and decoded per voxel into the five species by the
nife_intern pipeline (`results/fish_3d/ic/phi_xyz_<CH|DH>_d<day>.npy`, one field of view per condition and
day; shape 128 x 128 x n_z x 5, DH day 6 also at 256 x 256). The arrays are measured data of Heine et al.
and are read from that repository at run time; nothing of them is copied here. The repository holds only the
computed stress fields, the aggregated numbers below, the coarse section of the stress field in the figure
and, in the video, the phi = 0.3 surface of one stack at 2.9 um resolution.

- Values: each raw channel was divided by its own 99.5th percentile before the species decode, so the
  channel sum S is a dimensionless intensity (median 0.02 to 0.6, 99th percentile about 2, maximum up to
  18). It is not a calibrated volume fraction.
- phi: S is smoothed (Gaussian, 1.44 um laterally, 1 um in z, against the decode speckle) and mapped by
  phi = clip((S - S_10) / (S_95 - S_10), 0, 1): the 10th percentile of the stack is the background (every
  stack holds fluid or substratum), the 95th percentile is phi = 1. Film = phi > 0.3. This is an
  assumption; the sensitivity is in the checks.
- Voxel size: 0.18 um pixels (nife_intern notes), 1024 pixels per field of view, stride 8 in the saved
  arrays: 1.44 um laterally (0.72 um for the 256-wide stack), field 184 x 184 um. The z step is not stored
  with the arrays; it is taken as the pooled stack depth of that condition and day
  (`results/diffusion_fit/zprofiles_all_ti.csv`) divided by the number of slices of the stack, 0.8 to 3.5 um
  (nominal 2 um). The Leica metadata would settle it; the value scales the z axis, the film thickness H and
  with it beta.
- Substratum: in every stack the film has one sharp, flat edge and a ragged side with protrusions. The sharp
  edge is taken as the titanium disc, at the low z index (assumption: the stack starts at the substratum;
  the nife_intern analyses use both conventions). The base z_b(x, y) = first slice with phi > 0.3 in the
  column; a plane is fitted (least squares, residuals beyond two slices dropped). The stacks are tilted by
  0 to 15 degrees (DH day 10). Each column is shifted by an integer number of slices so that the plane is
  z = 0, slices below are dropped, and the box is cut at H_top + 30 % (H_top = top of the last slice where
  the film covers 5 % of the area). The mirror assumption is run as a check.
- Thickness: H = biovolume per substratum area (film voxels x voxel volume / area; the COMSTAT mean
  thickness), 5 to 31 um; local thickness of a column = its film voxels x dz.

| stack | dz (um) | tilt (deg) | H (um) | H_top (um) | film (% of box) | beta (mm^2/T*) |
|---|---|---|---|---|---|---|
| CH day 1 | 2.90 | 0.3 | 18.2 | 23.2 | 57 | 1.06e-04 |
| CH day 6 | 2.00 | 4.3 | 26.6 | 36.0 | 58 | 2.27e-04 |
| CH day 10 | 3.55 | 1.2 | 5.1 | 7.1 | 29 | 8.43e-06 |
| CH day 15 | 0.79 | 6.4 | 7.5 | 17.3 | 33 | 1.78e-05 |
| CH day 21 | 2.79 | 0.8 | 11.7 | 25.1 | 35 | 4.39e-05 |
| DH day 1 | 1.73 | 0.0 | 7.1 | 10.4 | 45 | 1.59e-05 |
| DH day 6 | 1.79 | 3.2 | 9.7 | 25.0 | 30 | 3.02e-05 |
| DH day 10 | 1.92 | 15.4 | 30.9 | 69.0 | 35 | 3.05e-04 |
| DH day 15 | 2.30 | 1.6 | 24.9 | 36.8 | 68 | 1.98e-04 |
| DH day 21 | 1.84 | 7.4 | 20.1 | 38.6 | 41 | 1.29e-04 |

## Model and scaling assumptions

- Growth as in `wp3_reference.py`: d phi/dt = beta lap phi + k_alpha alpha - P (max(0, phi - 1) + min(0, phi)),
  d alpha/dt = k_alpha phi, k_alpha = 1e-3 /T*, P = 5, T* = 1.1, RK4, zero flux on every face of the box
  (substratum, sides, top), phi(0) = the CLSM phi, alpha(0) = 1.
- beta = 0.02 mm^2/T* (H / 0.25 mm)^2: the spreading length relative to the film thickness is that of the
  thesis runs on the 0.25 mm film (sqrt(beta T*) = 0.6 H), the scaling the molar runs use. Assumption. With
  it the spreading length (2 to 20 um) is larger than the lateral structures of the stacks (a few voxels),
  so the field at T* is a smoothed version of the measurement (phi max falls from 1 to 0.5 to 0.9). The run
  with beta = 0 keeps the structure as measured and is the one that shows what the structure does.
- Mechanics at T*: small strain, hex8 with B-bar on the voxel mesh (elements 2.88 x 2.88 x dz um on the
  main mesh, lateral block mean of 2 x 2 voxels), E = 10 Pa (phi^2 + 1e-3) with phi at T*, nu = 0.49,
  eigenstrain (alpha - 1) I, substratum clamped, top free, lateral faces on rollers (normal displacement
  zero, the mirror symmetry of the zero-flux growth field; free and periodic faces in the checks). The
  material above and beside the film (phi below 0.3) stays in the mesh with its own small stiffness. pyamg
  smoothed aggregation + CG, relative residual 1e-9, 30 to 80 iterations, 0.06 to 0.57 M dof, about one
  minute per stack. The solver reproduces the uniform clamped layer (sigma = -E (alpha - 1)/(1 - nu)) to
  0.1 %.
- Uniform film: the same stack averaged in x and y, phi(z), grown and solved on a 2 x 2 column with rollers
  (the exact 1D solution); its band means are weighted with the real film's voxel count per slice.
- Measures over the film (phi(0) > 0.3): mean von Mises, bottom quarter / middle half / top quarter of
  H_top, 95th percentile, mean pressure p = -tr(sigma)/3 (compression positive), Spearman correlation of the
  local von Mises with the local phi and with the local thickness of the column. Maxima are listed but do
  not converge with the mesh (WALL_CONVERGENCE.md).

## Mesh check

Lateral 2.88 um (main mesh) against 1.44 um (every voxel), same z; and the central 92 x 92 um of the
256-wide stack of DH day 6 at 2.88, 1.44 and 0.72 um (the data at 0.72 um block-averaged, so the same
structure on three meshes). Von Mises in mPa (mean, bottom, middle, top, p95, p mean, max).

With spreading (beta as above), the field is smooth on every mesh:

| case | mean | bottom | middle | top | p95 | p mean | max |
|---|---|---|---|---|---|---|---|
| DH day 6, 2.88 um | 2.56 | 3.69 | 2.14 | 1.84 | 4.99 | 1.67 | 7.32 |
| DH day 6, 1.44 um | 2.56 | 3.69 | 2.14 | 1.84 | 5.00 | 1.67 | 7.34 |
| CH day 6, 2.88 um | 2.40 | 2.83 | 2.27 | 1.96 | 3.73 | 1.56 | 5.54 |
| CH day 6, 1.44 um | 2.40 | 2.83 | 2.27 | 1.96 | 3.72 | 1.56 | 5.54 |
| sub-region 2.88 um | 3.01 | 3.68 | 2.19 | 0.04 | 6.22 | 1.89 | 8.14 |
| sub-region 1.44 um | 3.01 | 3.69 | 2.20 | 0.05 | 6.26 | 1.89 | 8.19 |
| sub-region 0.72 um | 3.02 | 3.69 | 2.20 | 0.05 | 6.27 | 1.90 | 8.20 |

Without spreading (beta = 0), the structure is at the voxel scale and the means still move by about 1 %:

| case | mean | bottom | middle | top | p95 | p mean | max |
|---|---|---|---|---|---|---|---|
| DH day 6, 2.88 um | 8.68 | 10.53 | 7.96 | 7.70 | 19.89 | 5.13 | 25.40 |
| DH day 6, 1.44 um | 8.76 | 10.64 | 8.03 | 7.73 | 20.10 | 5.14 | 26.33 |
| CH day 6, 2.88 um | 7.05 | 2.63 | 9.25 | 6.66 | 16.95 | 4.35 | 23.29 |
| CH day 6, 1.44 um | 7.07 | 2.66 | 9.26 | 6.65 | 17.00 | 4.35 | 24.91 |
| sub-region 2.88 um | 10.49 | 11.34 | 9.48 | 0.19 | 21.42 | 6.19 | 25.49 |
| sub-region 1.44 um | 10.64 | 11.48 | 9.64 | 0.50 | 21.55 | 6.26 | 26.20 |
| sub-region 0.72 um | 10.70 | 11.54 | 9.71 | 0.52 | 21.61 | 6.28 | 27.11 |

The means and the 95th percentile change by 0.3 % or less with spreading and by 1 % (2.88 to 1.44 um) and
0.6 % (1.44 to 0.72 um) without, so the main mesh is sufficient for means and percentiles; the maximum grows
with refinement (3 to 6 % per step), as in the wall study.

## Results (main mesh, von Mises in mPa, real structure / uniform film with the same mean profile)

With spreading (beta scaled to H):

| stack | mean | bottom | middle | top | p95 | p mean | max | r_S(phi) | r_S(thickness) |
|---|---|---|---|---|---|---|---|---|---|
| CH day 1 | 2.69 / 2.64 | 3.71 / 3.67 | 2.41 / 2.37 | 1.15 / 0.91 | 4.20 / 3.76 | 1.79 / 1.76 | 4.78 / 3.76 | +0.36 | -0.14 |
| CH day 6 | 2.40 / 2.04 | 2.83 / 2.45 | 2.27 / 2.01 | 1.96 / 1.04 | 3.73 / 2.46 | 1.56 / 1.36 | 5.54 / 2.46 | +0.01 | +0.65 |
| CH day 10 | 3.52 / 2.78 | 4.42 / 3.67 | - / - | 2.08 / 1.35 | 7.98 / 3.67 | 2.34 / 1.85 | 12.83 / 3.67 | +0.78 | -0.44 |
| CH day 15 | 3.30 / 2.12 | 4.65 / 3.01 | 2.67 / 1.76 | 1.23 / 0.21 | 6.92 / 3.13 | 2.17 / 1.42 | 10.21 / 3.14 | +0.55 | -0.13 |
| CH day 21 | 2.48 / 2.05 | 3.69 / 3.16 | 1.79 / 1.42 | 0.43 / 0.06 | 5.00 / 3.39 | 1.65 / 1.37 | 7.74 / 3.39 | +0.44 | -0.03 |
| DH day 1 | 3.21 / 2.67 | 4.34 / 3.80 | 3.14 / 2.69 | 1.67 / 0.70 | 6.01 / 3.80 | 2.11 / 1.78 | 12.82 / 3.80 | +0.41 | +0.37 |
| DH day 6 | 2.56 / 1.17 | 3.69 / 1.69 | 2.14 / 1.09 | 1.84 / 0.22 | 4.99 / 1.73 | 1.67 / 0.78 | 7.32 / 1.73 | +0.43 | -0.19 |
| DH day 10 | 1.69 / 1.53 | 2.39 / 2.16 | 1.24 / 1.16 | 0.16 / 0.08 | 3.32 / 2.27 | 1.12 / 1.02 | 3.83 / 2.28 | +0.53 | +0.05 |
| DH day 15 | 2.74 / 2.49 | 3.21 / 2.84 | 2.59 / 2.43 | 2.28 / 1.57 | 4.36 / 2.85 | 1.81 / 1.66 | 5.37 / 2.85 | +0.45 | +0.40 |
| DH day 21 | 2.06 / 1.59 | 2.57 / 1.83 | 1.99 / 1.61 | 1.21 / 0.72 | 3.48 / 1.86 | 1.37 / 1.06 | 5.22 / 1.86 | +0.18 | +0.34 |

Without spreading (beta = 0, the structure as measured):

| stack | mean | bottom | middle | top | p95 | p mean | max | r_S(phi) | r_S(thickness) |
|---|---|---|---|---|---|---|---|---|---|
| CH day 1 | 7.39 / 6.12 | 11.82 / 10.33 | 5.86 / 4.84 | 2.52 / 0.39 | 20.25 / 15.36 | 4.79 / 4.08 | 26.70 / 15.36 | +0.95 | -0.24 |
| CH day 6 | 7.05 / 5.07 | 2.63 / 1.78 | 9.25 / 7.36 | 6.66 / 0.88 | 16.95 / 10.42 | 4.35 / 3.38 | 23.29 / 10.54 | +0.94 | +0.02 |
| CH day 10 | 9.89 / 4.74 | 10.95 / 6.39 | - / - | 8.20 / 2.10 | 21.35 / 6.39 | 6.50 / 3.16 | 22.70 / 6.39 | +0.96 | -0.50 |
| CH day 15 | 9.24 / 3.29 | 11.28 / 4.13 | 8.57 / 3.08 | 3.10 / 0.09 | 20.67 / 5.45 | 5.79 / 2.19 | 25.16 / 5.46 | +0.91 | -0.25 |
| CH day 21 | 7.95 / 4.03 | 12.82 / 7.61 | 4.86 / 1.77 | 2.76 / 0.02 | 21.30 / 9.24 | 5.09 / 2.69 | 24.76 / 9.24 | +0.94 | -0.23 |
| DH day 1 | 8.51 / 6.48 | 4.84 / 3.36 | 10.89 / 8.92 | 3.64 / 0.39 | 19.80 / 15.57 | 5.47 / 4.32 | 26.86 / 15.57 | +0.95 | -0.03 |
| DH day 6 | 8.68 / 1.48 | 10.53 / 2.02 | 7.96 / 1.46 | 7.70 / 0.16 | 19.89 / 2.28 | 5.13 / 0.98 | 25.40 / 2.28 | +0.91 | -0.23 |
| DH day 10 | 5.39 / 3.09 | 6.66 / 4.01 | 4.84 / 2.64 | 0.54 / 0.03 | 17.00 / 6.23 | 3.22 / 2.06 | 26.37 / 6.26 | +0.92 | -0.10 |
| DH day 15 | 6.15 / 4.93 | 3.91 / 2.56 | 7.59 / 6.46 | 1.82 / 0.35 | 15.66 / 9.85 | 3.85 / 3.28 | 23.23 / 9.90 | +0.94 | +0.06 |
| DH day 21 | 6.67 / 2.66 | 5.38 / 1.25 | 7.53 / 3.45 | 3.60 / 0.46 | 15.67 / 4.39 | 3.88 / 1.77 | 24.71 / 4.41 | +0.90 | -0.08 |

CH day 10 has only two film slices (dz = 3.55 um, H_top = 7.1 um), so its middle band is empty. The
reference stress of a fully constrained uniform layer with phi = 1 is E (alpha - 1)/(1 - nu) = 21.6 mPa.

## What the real structure changes against a uniform film

The real structure raises the mean von Mises of the film above the uniform film with the same profile in
every stack: by 2 to 119 % (median 21 %) with the scaled spreading and by 21 to 488 % (median 86 %) without
it, and the 95th percentile by 12 to 188 % and 27 to 774 %. The reason is the stiffness weighting E ~ phi^2:
the uniform film carries the growth strain in a soft, diluted layer, while the real film concentrates both
the stiffness and the eigenstrain in dense voxels, and the local von Mises follows the local phi (Spearman
+0.9 without spreading, +0.2 to +0.8 with it). The local thickness of the column matters much less (|r_S|
mostly below 0.3, sign varying), so the stress is set by the local density, not by how thick the film is
above a point. Dense clumps near the top of the film carry the largest gain against the uniform film (top
quarter 4 to 170 times), the band on the substratum the smallest (1.1 to 5 times), because the uniform film
is already most constrained there. The sparse DH stacks (days 6 and 21, film 30 to 40 % of the box) gain the
most and the dense ones (CH days 1 and 6, DH days 1 and 15) the least. With the scaled spreading the field
at T* is smooth and the difference shrinks to tens of per cent; the structure effect is a property of the
measured film, which the growth model smooths out within one T*. Between the conditions, without spreading
the CH stacks carry 8.3 mPa on average and the DH stacks 7.1 mPa (7.0 to 9.9 against 5.4 to 8.7); with
spreading 2.9 against 2.5 mPa; there is no trend with the day, the stack-to-stack scatter (one field of view
per day) is larger than any difference between the conditions.

## Checks (DH day 6, main mesh, with spreading unless stated)

- Lateral faces: rollers 2.56 mPa film mean, periodic 2.46, free 1.98; over the interior (columns in the
  central half of the width) 2.96 / 2.92 / 2.78 mPa, so free faces lower the interior mean by 6 % (beta = 0:
  10.15 / 9.84 / 8.66, 15 %) and periodic faces by 1 to 3 %. The rollers are kept.
- phi map: (S_10, S_95) 2.56 mPa; (0, S_95) 2.41; (S_10, S_99) 1.86; (S_10, S_90) 2.89. The map sets the
  level of phi and so the stress by up to 30 %; the ratio to the uniform film and the correlations change
  little (r_S(phi) +0.3 to +0.5). Without the background subtraction the noise voxels in the fluid pass the
  threshold and H_top rises from 25 to 46 um.
- Mirror orientation (substratum at the high z index): 1.94 mPa against 2.56 (uniform 0.94 against 1.17);
  the ratio to the uniform film stays at about 2. The orientation changes the level, not the conclusion.

## Limitations

- One field of view per condition and day; the stacks differ in z step, depth and background level, so
  the differences between days and conditions are within the scatter of the fields of view.
- The intensity-to-phi map and the smoothing are assumptions; the stress level scales with them (up to
  30 %), the ratios and correlations do not.
- The z step is reconstructed from the pooled depth of the day, 0.8 to 3.5 um, and the orientation of the
  stack (substratum at the low index) is inferred from the shape of the film; both need the Leica metadata.
- The substratum is a fitted plane with integer-slice shifts (one to two slices of error); DH day 10 is
  tilted by 15 degrees and its upper 40 um are sparse.
- Growth to T* = 1.1 only, small strain, E = 10 Pa; the stresses are in mPa and scale linearly with E and
  with the growth. The scaled spreading smooths the measured structure within one T*; the beta = 0 run is
  the structure as measured.
- The maximum does not converge with the mesh (voxel corners); the means, bands and percentiles do (1 %).
