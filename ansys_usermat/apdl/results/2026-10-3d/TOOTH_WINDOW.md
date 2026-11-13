# Cervical window on the molar crown at physical size (7 Nov 2026)

Set-up: `papers/coupling_paper/REAL_TOOTH_WINDOW.ja.md`. The code works in the units of the 2 mm cube (Klempt 2024's
20 um cube x 100), so 1 code mm = 10 um and the constants of Klempt 2024 Table 2 are used unchanged. The crown STL
(`tier2b_real/_crown_real_watertight.stl`) is scaled by 100 and the window is a cube of edge W = 30 code mm = 300 um.

## Commands

Fields and point model (the standard options of the paper runs plus the window options; `ECO_THREADS=4`):

    python ansys_usermat/apdl/coupled_living_consumption.py --nf 64 --front 1 --phi-min 0.001 --v2 --fortran \
        --inherit-living --alpha0 100 --gab 1 --gc 1 --tab 0.333333 --tvd --limiter superbee --transient-c 1e8 --amg \
        --obstacle stl --stl-scale 100 --box 30 --window mixed --only window_ab
    (--window layered for the layered start; --only window for the run without antibiotic; --nf 32 / 128)

New options (`coupled_living_consumption.py`): `--window mixed|layered` (the example "window", keys `window`,
`window_ab`), `--stl-scale S` (the STL scaled by S), `--film-h0 H` (start film thickness in code mm, default 2 = 20 um).
Tag: `..._obststl_win<comp>_stlx100_box30_amg_window_ab`; outputs `coupled_living_consumption_<tag>.json` here and
`figdata/coupled_living_window_ab_<tag>.npz` (snapshots at T* = 0.35, the first outer step after the onset, and 1).

Film stress (rigid crown, roller faces; `real_tooth_submodel.solve_submodel(..., roller=...)` is new, the default
path is unchanged, checked against the committed version on the 40^3 crown run: max |dvM| = 3e-13 Pa):

    python ansys_usermat/apdl/real_tooth_coupled.py --window --nf 64 --box 30 --cases window_ab \
        --opts front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_superbee_wgmirror_obststl_winmixed_stlx100_box30_amg
    -> real_tooth_coupled_transc1e+08_tvd_superbee_wgmirror_obststl_win<comp>_stlx100_box30_amg.json (keys window_ab_nf<NF>),
       figdata/real_tooth_coupled_stress_window_ab_nf<NF>_T<T>_h<1|2>.npz

Figures: `fig_tooth_window_setup.py` (set-up, `assets/fig_tooth_window_setup.pdf/.png`) and
`fig_tooth_window_section.py --nf 64 --window mixed [--case window]` (section y = 0 at the start and at T* = 1,
`assets/fig_tooth_window_section_<comp>_nf<NF>[_noab].png`; drawn for the three 64^3 runs and the 128^3 run).

## Geometry and boundary conditions as implemented

- Window placement: the crown's lowest z (the flat cap of the cervical base, z = 32.5 mm in the STL) is the bottom
  face z = -W/2; the window centre is the point of the crown surface with the smallest x at W/2 above the base
  (STL point (-72.139, -40.421, 32.65) mm, scaled x 100). The crown occupies the +x side: 48.1 % of the cells
  (126 039 of 262 144 at 64^3); at y = 0 the surface is at x = 0 in the middle, x = +33 um at the bottom and x = -19 um
  at the top row (the flank leans outward going up, dx/dz = -0.154 between z = -W/4 and +W/4, i.e. 81 degrees from the
  horizontal: roughly vertical). Across y the surface moves by about 4 um (the crown's radius of curvature there is
  about 3 mm). The crown reaches the bottom face (39 % of the bottom cell layer is tooth, 56 % of the top layer).
- SDF sign: pyvista's implicit distance gives the wrong sign where the cap is the closest surface (a layer of
  "outside" cells inside the crown just above the cap); the window example takes the sign from the point-in-closed-
  surface test. **The crown example (`--obstacle stl --box 10`, `tooth_ab`) still uses the implicit sign: at 40^3,
  1397 of 8596 inside cells (z indices 0-7, up to 2 mm above the base) are classified as open there.** Not changed
  here (existing results); to be looked at.
- Start film: cells with 0 <= SDF < 2 code mm: 17 911 cells at 64^3, 4.27 cells across (2.13 at 32^3, 8.5 at 128^3),
  mean thickness 20.5 um by column count; mixed start with Fritsch's shares (0.25, 0.30), layered: A on SDF < 1,
  B on 1 <= SDF < 2.
- Field BCs (checked on the matrix of `Reaction(q, "xminus")` at 16^3): the Dirichlet ghost term (2/h^2) and the
  right-hand side only in the cell layer x = -W/2, no boundary term in any other face layer (zero flux), tooth rows
  u = 0; with q = 0 the solution is c = 1 to 7e-13. In the runs, c in the first cell layer at x = -W/2 is 0.982-0.988
  (face value 1 by the ghost cell; the cell centre is h/2 inside), a there 0.981 after the onset; c = 0 in every tooth
  cell. The antibiotic from T* = 1/3 on the same face.
- Mechanics BCs: every node of a tooth element u = 0; u_y = 0 on the mesh faces at y = +-W/2, u_z = 0 at z = +-W/2
  (1598 / 2469 roller nodes at T* = 0.35 / 1 at 64^3); the face x = -W/2 and the film surface free (elements with
  phi <= 0.05 dropped). The crop box reaches all four roller faces (the film covers the whole crown surface in the
  window).
- phi_3D stays in [0, 1.0002] over every run (penalty overshoot 2e-4, as the cube examples).

## Numbers (means over the film, phi_3D > 0.5 and not tooth; stress h = h_f, in brackets h = 2 h_f)

| run | T* | film cells | share A | sum phibar | c (min) | a | alpha - 1 | thickness mean (min-max) [um] | vM film [mPa] | vM p95 / max [mPa] |
|---|---|---|---|---|---|---|---|---|---|---|
| 32^3 mixed, antibiotic | 0.35 | 2 233 | 0.453 | 0.496 | 0.046 (0.012) | 0.045 | 3.0e-4 | 27 (sdf max 29) | - | - |
| 32^3 mixed, antibiotic | 1 | 4 319 | 0.421 | 0.037 | 0.340 (0.247) | 0.331 | 7.0e-4 | 40 (28-56) | - | - |
| 64^3 mixed, antibiotic | 0.35 | 21 816 | 0.453 | 0.498 | 0.055 (0.013) | 0.054 | 3.1e-4 | 25.0 (18.8-28.1) | 5.92 (4.90) | 8.55 / 12.7 |
| 64^3 mixed, antibiotic | 1 | 34 157 | 0.392 | 0.049 | 0.247 (0.184) | 0.232 | 7.0e-4 | 39.1 (32.8-46.9) | 13.87 (12.35) | 23.5 / 35.7 |
| 64^3 layered, antibiotic | 0.35 | 21 768 | 0.419 | 0.538 | 0.053 (0.012) | 0.052 | 3.1e-4 | 24.9 (18.8-28.1) | 5.92 (4.90) | 8.56 / 12.7 |
| 64^3 layered, antibiotic | 1 | 27 146 | 0.344 | 0.295 | 0.076 (0.017) | 0.074 | 8.1e-4 | 31.1 (23.4-37.5) | 15.56 (13.52) | 23.9 / 36.3 |
| 64^3 mixed, no antibiotic | 0.35 | 21 816 | 0.455 | 0.536 | 0.055 (0.013) | 0 | 3.1e-4 | 25.0 (18.8-28.1) | 5.92 (4.90) | 8.55 / 12.7 |
| 64^3 mixed, no antibiotic | 1 | 26 027 | 0.466 | 0.580 | 0.039 (0.008) | 0 | 8.3e-4 | 29.8 (23.4-37.5) | 15.89 (13.73) | 24.0 / 36.3 |
| 128^3 mixed, antibiotic | 0.35 | 166 554 | 0.453 | 0.497 | 0.054 (0.013) | 0.054 | 3.2e-4 | 23.8 (21.1-25.8) | 6.26 (5.77) | 8.14 / 12.3 |
| 128^3 mixed, antibiotic | 1 | 256 998 | 0.358 | 0.102 | 0.160 (0.118) | 0.155 | 7.3e-4 | 36.8 (30.5-42.2) | 14.46 (13.67) | 22.5 / 34.5 |

Thickness = number of film cells along x per (y, z) column times h (the flank is nearly vertical, so the column
count is the thickness to within 2 %); the volume over the exposed voxel area of the crown gives 20-25 % less
(the staircase area overcounts). Von Mises by distance from the crown at T* = 1, 64^3 mixed with antibiotic:
0-10 um 21.2, 10-20 um 19.3, 20-30 um 10.6, 30-40 um 2.3 mPa (the old layers carry alpha - 1, the new ones hardly
any). The stress maximum sits at the bottom corner of the window (x = +26, y = -134, z = -134 um), where the roller
faces y = -W/2 and z = -W/2 meet the crown: a corner of the model, not of the film.

At the antibiotic onset the three 64^3 runs agree (same start, the antibiotic has not acted yet); the mixed and the
layered start differ only in the composition. Reading of the T* = 1 numbers:

- With antibiotic the living fraction of the mixed film falls from 0.50 to 0.05; the layered film keeps 0.30 (the
  inner A layer stays living: the antibiotic is consumed in the outer layer, a in the film 0.07 against 0.23).
- The film with antibiotic grows thicker (39 um) than without (30 um): the dead film consumes less nutrient, c in the
  film rises from 0.04 to 0.25, and the front speed r c/(k + c) with it (the thesis limitation that dead biofilm
  keeps spreading).
- The film does not reach the window faces at T* = 1 (thickness at most 47 um in a 300 um window); the film is thin
  compared with the window, so the cut-face and gingiva conditions act only through the stress of the film ends.

## Timing and memory (this cloud container, 4 cores, ECO_THREADS=4)

| grid | run | wall | max RSS |
|---|---|---|---|
| 32^3 | mixed, antibiotic | 15 s | 0.43 GB |
| 64^3 | mixed, antibiotic | 120 s | 0.76 GB |
| 64^3 | layered, antibiotic | 101 s | 0.76 GB |
| 64^3 | mixed, no antibiotic | 71 s | 0.69 GB |
| 128^3 | mixed, antibiotic | 1460 s (24 min) | 3.2 GB |

Stress at 64^3: 6-10 s per solve (79k-122k dof, 44-54 CG iterations), 1 s at 2 h_f; at 128^3 55-87 s
(569k-871k dof, 5843 / 8887 roller nodes), 7-10 s at 2 h_f. The 128^3 checkpoint (every 120 s) is 290 MB.
192^3 would take about 3.5 h and 10 GB by the same scaling.

## What is not converged / open

- Grid (mixed with antibiotic, T* = 1): thickness 40 / 39 / 37 um, film c 0.34 / 0.25 / 0.16, living fraction
  0.037 / 0.049 / 0.102, film vM (h = h_f) - / 13.9 / 14.5 mPa at 32^3 / 64^3 / 128^3. The thickness and the stress
  are within 5 % between 64^3 and 128^3; the film's c and living fraction are not converged (the nutrient depletion
  length is 1 code mm = 10 um = 2 cells at 64^3, 4 cells at 128^3, and the living fraction depends on c through the
  point model), so the composition numbers of the window need 128^3 or finer. At the onset (T* = 0.35) all three
  grids agree to 2 %. 192^3 (as the spec) is not run here (about 3.5 h).
- The von Mises by distance at 128^3 (0-10 / 10-20 / 20-30 / 30-40 um: 20.4 / 19.6 / 10.6 / 2.1 mPa) matches 64^3
  (21.2 / 19.3 / 10.6 / 2.3) to within 4 %.
- The stress at h = h_f vs 2 h_f differs by 11-20 % (the first element layer on the voxel staircase; the same as the
  crown example). The maximum is at a model corner.
- The front speed: the spec expects about 100 um/T* (c = 1). Here c at the film surface is about 0.5 and in the film
  0.05-0.25, because the nutrient is held at c = 1 only on the face 280 um away and diffuses through the saliva region
  without flow (d = 1e8 code mm^2/T*, consumption g = 1e8/T* x phi sum phibar: the drop across the 280 um of saliva is
  of the same order as the film's uptake), so the film grows by 10-20 um in T* = 1, not 100 um. If the saliva is meant
  to be well mixed, the nutrient should be held closer to the film (e.g. a smaller window in x, or c = 1 on the film
  surface); a decision for the paper.
- The nutrient depletion length (10 um) is resolved by 2 cells at 64^3; the min c in the film at the crown (0.01 at
  T* = 0.35) is a few cells from the surface: 128^3 (4 cells) is the first grid that resolves it.

## Caption of the set-up figure (`assets/fig_tooth_window_setup.pdf`)

Cervical window on the molar crown at physical size (1 code mm = 10 um). (a) The crown (Open-Full-Jaw,
Gholamalizadeh et al. 2022) with the window (W = 300 um) at the -x flank of the cervical region; the bottom face of
the window lies in the plane of the cervical base. (b) The window with the crown as the voxel body of the 64^3 grid
(grey), the start film of 20 um (mixed start) and the boundary conditions: c = 1 and, from T* = 1/3, a = 1 on the
saliva face x = -W/2 (free in the mechanics); zero flux and roller (u_n = 0) on the cut faces y = +-W/2 and z = +W/2;
zero flux and u_z = 0 on the gingiva face z = -W/2; the crown rigid (u = 0), no biofilm and zero flux. The octant
x < 0, y < 0, z > 0 is removed; z is vertical (the camera of the cube examples with y and z exchanged, so that the open
side faces the viewer). (c) The layered start: species A on the inner 10 um, B on the outer 10 um. (d) The section
y = 0 at T* = 0 with the boundary conditions on its edges; scale bar 100 um.

---

# Variant 2 (7 Nov 2026, the user's decision): the saliva well mixed, `--saliva-mixed 0.05`

The results above hold the nutrient only on the face x = -W/2 (variant 1, diffusion through 280 um of still
saliva; kept for comparison). The paper uses this variant: the saliva is well mixed, so at the start of every outer
step the open cells (not tooth) with phi_3D < phi_s = 0.05 hold c = 1 and, from T* = 1/3, a = 1 (Dirichlet cells, the
fixed-cell mechanism of `--ansys-bc`), in addition to the face x = -W/2. A boundary condition that follows the film
surface; the equations are unchanged. Tag `_salmix0.05` (`hist["saliva_mixed"]`, `hist["saliva_cells"]` per step).
The AMG hierarchy of the field solves is rebuilt whenever the fixed set changes (every outer step here).

    python ansys_usermat/apdl/coupled_living_consumption.py --nf 64 --front 1 --phi-min 0.001 --v2 --fortran \
        --inherit-living --alpha0 100 --gab 1 --gc 1 --tab 0.333333 --tvd --limiter superbee --transient-c 1e8 --amg \
        --obstacle stl --stl-scale 100 --box 30 --saliva-mixed 0.05 --window mixed --only window_ab
    python ansys_usermat/apdl/real_tooth_coupled.py --window --nf 64 --box 30 --cases window_ab \
        --opts front1_phimin0.001_v2_fortran_inhliv_a100_gab1_gc1_tab0.333333_transc1e+08_tvd_superbee_wgmirror_obststl_winmixed_stlx100_salmix0.05_box30_amg
    python ansys_usermat/apdl/fig_tooth_window_section.py --nf 64 --window mixed --saliva-mixed 0.05 [--case window]
        -> assets/fig_tooth_window_section_<comp>_nf<NF>_salmix0.05[_noab].png

**phi_s sensitivity (64^3, mixed with antibiotic, T* = 1):** phi_s = 0.05: thickness 61.4 um, living 0.201, c 0.224,
share A 0.606, film vM 13.11 mPa; phi_s = 0.01: 60.0 um, 0.191, 0.221, 0.594, 13.21 mPa, i.e. 2 % in the thickness and
5 % in the living fraction (the cells between phi = 0.01 and 0.05 are the outermost layer of the front, where the
nutrient is already 0.85). phi_s = 0.05 is used (the same threshold as the dropped elements of the stress mesh).

**Does the film feel the saliva face?** No. At T* = 1 the film (phi_3D > 0.5) is 61 um thick (max 70 um) and the
colonised front (phi_3D > 0.05) reaches 74 um from the crown; the saliva face is 150-180 um from the crown. The front
is diffuse (phi decays from 0.5 to 0.05 over about 15-25 um); the stress mesh (phi > 0.05) extends to the front. No
W = 400 um run was needed. The front moves at 30-40 um/T* (not 100 um/T*): c at the phi = 0.5 cells is 0.5-0.6 (0.76 in
the front cells, 0.2 in the film), so r c/(k + c) is 3-4 mm/T* there, and the thickness is measured at phi = 0.5.

## Numbers, variant 2 (film = phi_3D > 0.5 and not tooth; thickness by column count; stress h = h_f, [2 h_f])

| run | T* | film cells | share A | sum phibar | c (min) | a | alpha - 1 | thickness um (min-max) | vM film [mPa] | vM p95 / max |
|---|---|---|---|---|---|---|---|---|---|---|
| 64^3 mixed, antibiotic | 0.35 | 30 963 | 0.492 | 0.653 | 0.214 (0.037) | 0.208 | 2.8e-4 | 35.4 (32.8-42.2) | 5.34 [4.72] | 8.27 / 12.5 |
| 64^3 mixed, antibiotic | 1 | 53 676 | 0.606 | 0.201 | 0.224 (0.120) | 0.217 | 6.7e-4 | 61.4 (56.2-70.3) | 13.11 [12.14] | 22.7 / 34.4 |
| 64^3 layered, antibiotic | 0.35 | 30 958 | 0.301 | 0.714 | 0.211 (0.034) | 0.205 | 2.8e-4 | 35.4 (32.8-42.2) | 5.33 [4.72] | 8.27 / 12.5 |
| 64^3 layered, antibiotic | 1 | 52 125 | 0.195 | 0.307 | 0.198 (0.013) | 0.195 | 6.7e-4 | 59.7 (46.9-70.3) | 12.69 [11.70] | 22.8 / 34.7 |
| 64^3 mixed, no antibiotic | 0.35 | 30 963 | 0.503 | 0.741 | 0.214 (0.037) | 0 | 2.8e-4 | 35.4 (32.8-42.2) | 5.34 [4.72] | 8.27 / 12.5 |
| 64^3 mixed, no antibiotic | 1 | 47 881 | 0.556 | 0.808 | 0.062 (0.002) | 0 | 7.0e-4 | 54.8 (46.9-65.6) | 13.06 [11.98] | 22.9 / 34.9 |
| 128^3 mixed, antibiotic | 0.35 | 257 930 | 0.497 | 0.677 | 0.237 (0.036) | 0.234 | 2.7e-4 | 36.9 (32.8-42.2) | 5.32 [5.04] | 7.85 / 12.0 |
| 128^3 mixed, antibiotic | 1 | 449 476 | 0.647 | 0.284 | 0.195 (0.085) | 0.192 | 6.6e-4 | 64.3 (58.6-75.0) | 13.03 [12.57] | 21.6 / 34.0 |
| 128^3 layered, antibiotic | 0.35 | 257 607 | 0.288 | 0.745 | 0.234 (0.033) | 0.231 | 2.7e-4 | 36.9 (32.8-42.2) | 5.32 [5.04] | 7.85 / 12.0 |
| 128^3 layered, antibiotic | 1 | 446 978 | 0.183 | 0.321 | 0.224 (0.011) | 0.223 | 6.6e-4 | 63.9 (53.9-75.0) | 12.81 [12.35] | 21.6 / 34.0 |
| 128^3 mixed, no antibiotic | 0.35 | 257 930 | 0.509 | 0.769 | 0.237 (0.036) | 0 | 2.7e-04 | 36.9 (32.8-42.2) | 5.32 [5.04] | 7.9 / 12.0 |
| 128^3 mixed, no antibiotic | 1 | 422 140 | 0.583 | 0.827 | 0.087 (0.002) | 0 | 6.8e-04 | 60.4 (49.2-70.3) | 13.03 [12.55] | 21.7 / 34.0 |

Von Mises by distance from the crown at T* = 1 (mixed with antibiotic): 64^3 21.3 / 19.6 / 17.0 / 12.1 / 6.8 / 2.0 mPa,
128^3 20.5 / 19.8 / 17.4 / 12.9 / 8.2 / 2.8 mPa for 0-10 / 10-20 / ... / 50-60 um: the old layers (first 30 um) carry
alpha - 1 and 17-21 mPa, the new ones little. The maximum sits at the bottom corner of the window where the roller faces
y = -W/2 and z = -W/2 meet the crown (a model corner), as in variant 1.

**64^3 vs 128^3 (T* = 1):** thickness 61.4 -> 64.3 um (mixed, +5 %), 59.7 -> 63.9 (layered), 54.8 -> 60 (no
antibiotic); film vM 13.11 -> 13.03 mPa (-1 %), 12.69 -> 12.81 (layered); living fraction 0.20 -> 0.28 (mixed with
antibiotic), 0.31 -> 0.32 (layered), 0.81 -> 0.83 (no antibiotic); film c 0.224 -> 0.195 (mixed), 0.198 -> 0.224
(layered), 0.062 -> 0.087 (no antibiotic); share A 0.61 -> 0.65 (mixed), 0.195 -> 0.183 (layered). At the onset
(T* = 0.35) the two grids agree within 4 % in every column. The thickness and the stress are converged to about 5 %;
the living fraction of the mixed film with antibiotic (0.20 -> 0.28) is the least converged number (it depends on the
nutrient in the inner film, depletion length 10 um = 2 cells at 64^3, 4 at 128^3). 192^3 is not run (about 1.3 h per
run and 7 GB by the scaling of the table below).

Readings (both grids): the layered start keeps the inner A layer living (0.32 at T* = 1 against 0.28 mixed; the inner
layer's c_min 0.01 against 0.09, the antibiotic reaches it later), the share of A in the film grows with antibiotic in
the mixed start (0.50 -> 0.65: B dies faster) and falls in the layered start (the outer B layer is colonised first);
without antibiotic the film stays 0.83 living, consumes the nutrient (c 0.09 in the film) and grows more slowly (60
against 64 um), the known limitation that the dead film keeps spreading.

**Film thickness over time, 128^3** (as the 64^3 table above):

| run | T* = 0.1 | 0.25 | 0.35 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| 128^3 mixed, antibiotic | 25 (0.60) | 33 (0.72) | 37 (0.68) | 43 (0.12) | 55 (0.19) | 64 (0.28) |
| 128^3 layered, antibiotic | 25 (0.63) | 33 (0.74) | 37 (0.74) | 43 (0.59) | 53 (0.42) | 64 (0.32) |
| 128^3 mixed, no antibiotic | 25 (0.60) | 33 (0.72) | 37 (0.77) | 43 (0.80) | 52 (0.82) | 60 (0.83) |

## Timing and memory, variant 2 (4 cores, ECO_THREADS=4; the AMG hierarchy rebuilt every outer step)

| grid | run | wall | max RSS |
|---|---|---|---|
| 64^3 | mixed, antibiotic | 114 s | 0.67 GB |
| 64^3 | mixed, antibiotic, phi_s = 0.01 | 98 s | 0.67 GB |
| 64^3 | layered, antibiotic | 101 s | 0.67 GB |
| 64^3 | mixed, no antibiotic | 83 s | 0.64 GB |
| 128^3 | mixed, antibiotic | 1211 s (20 min) | 2.7 GB |
| 128^3 | layered, antibiotic | 1153 s | 2.7 GB |
| 128^3 | mixed, no antibiotic | 1079 s | 2.4 GB |

Stress (h = h_f): 64^3 10-16 s (116k-198k dof), 128^3 83-175 s (0.86-1.54 M dof, 50 CG iterations); 2 h_f 1-17 s.

Figures of variant 2: `assets/fig_tooth_window_section_{mixed,layered}_nf{64,128}_salmix0.05.png`,
`assets/fig_tooth_window_section_mixed_nf{64,128}_salmix0.05_noab.png`; the set-up figure
`assets/fig_tooth_window_setup.pdf/.png` now names the saliva condition (phi_3D < 0.05 and the face x = -W/2: c = 1,
a = 1 from T* = 1/3; free). Caption (b): "... c = 1 and, from T* = 1/3, a = 1 in the saliva (the open cells with
phi_3D < 0.05) and on the face x = -W/2 (free in the mechanics); ..." (the rest as the caption above).

**Film thickness over time** (column count, phi_3D > 0.5, in um; in brackets the film's living fraction sum_i phibar_i;
start 20 um; the antibiotic from T* = 1/3). For the later comparison with plaque: 15 um at 1 day, 50 um at 7 days
(Al-Ahmad et al. 2007), so T* = 1 of this run corresponds to a few days of growth from a 1-day film.

| run | T* = 0.1 | 0.25 | 0.35 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| 64^3 mixed, antibiotic | 25 (0.59) | 31 (0.68) | 35 (0.65) | 40 (0.12) | 53 (0.11) | 61 (0.20) |
| 64^3 layered, antibiotic | 25 (0.62) | 31 (0.71) | 35 (0.71) | 40 (0.58) | 50 (0.41) | 60 (0.31) |
| 64^3 mixed, no antibiotic | 25 (0.59) | 31 (0.68) | 35 (0.74) | 40 (0.78) | 49 (0.80) | 55 (0.81) |

## Variant 2 on 192^3 (mixed, antibiotic; cloud, 7 Nov)

Field 4879 s (83 min, 4 cores). Stress at h = h_f ran out of memory at T* = 1 (about 5 M dof); at 2 h_f it is fine.

| grid | thickness um | sum phibar | share A | alpha-1 [1e-4] | vM film mPa (2 h_f) |
|---|---|---|---|---|---|
| 64^3 | 61.4 | 0.201 | 0.606 | 6.70 | 12.14 |
| 128^3 | 64.3 | 0.284 | 0.647 | 6.60 | 12.57 |
| 192^3 | 65.8 | 0.320 | 0.666 | 6.58 | 12.67 |

64/128/192^3: vM p 1.55, extrapolated 12.79 mPa, GCI 1.2 %; alpha - 1 GCI 0.5 %; sum phibar p 0.54, extrapolated 0.47,
GCI 58 % (not converged: the 10 um depletion length is 2-6 cells); share A p 0.42, GCI 19 %; thickness by column count not in
the asymptotic range. 256^3 asked of fifa (CLOUD_TO_FIFA part 11, rewritten).

## Time scale from the film thickness (8 Nov)

`ansys_usermat/apdl/fig_tooth_window_thickness.py` -> `assets/fig_tooth_window_thickness.pdf/.png`. Thickness = film
(phi_3D > 0.5) volume per tooth area, from hist "film" "thickness_mean_mm" (code mm x 10 = um); it starts at 18 um
(the 20 um start film on the voxel wall) and grows nearly linearly. In situ plaque on enamel: about 15 um after one day and
50 um after seven (Al-Ahmad et al. 2007), i.e. about 5.8 um per day.

| run (variant 2) | slope T* = 0.5-1 [um/T*] | days per T* | T* at 50 um | thickness at T* = 1 [um] |
|---|---|---|---|---|
| 128^3 mixed, antibiotic | 32.5 | 5.6 | 0.96 | 50.5 |
| 192^3 mixed, antibiotic | 34.3 | 5.9 | 0.92 | 51.7 |
| 128^3 layered, antibiotic | 32.8 | 5.6 | 0.97 | 50.2 |
| 128^3 mixed, no antibiotic | 28.0 | 4.8 | 0.98 | 47.4 |

Reading: with the constants of Klempt et al. 2024 at their own size, the film grows like an in situ plaque if one T* is
about five to six days; the film of 18 um at the start then corresponds to about day 1.5 and the end (50 um) to about
day 7. Two straight lines through two measured points: an order of magnitude, not a calibration. (The earlier column-count
thickness of 64 um is the largest extent across the film, the volume per area is the mean.)

## Calibrated five species in the window (8 Nov)

`--eco-file abaqus_composition/eco_cal5_CH.txt` (the CH MAP of Chapter 3, Heine et al. 2025 data: c* = 25, alpha* = 0, so no
antibiotic), the film starting from the calibration's state at day 1 in every film cell (new branch in
coupled_living_consumption.py: window with --eco-file and the mixed composition), time link of the calibration
(T* = 1 = day 1 to day 21, --eco-steps 2262). `--v2` is left out (it only sets the cube's antibiotic and alpha*_0 = 50,
which the file's alpha* = 0 forbids). Results in `results/2026-10-validation/` (the solver writes --eco-file runs there).

```
ECO_THREADS=4 python3 coupled_living_consumption.py --nf 64 --front 1 --phi-min 0.001 --fortran --inherit-living --gc 1 --tvd \
  --limiter superbee --transient-c 1e8 --amg --obstacle stl --stl-scale 100 --box 30 --saliva-mixed 0.05 --window mixed \
  --eco-file ../../abaqus_composition/eco_cal5_CH.txt --only window          # 64^3: 2 min
```

Film living fractions phibar_i (So, An, Vei, Fn, Pg) and their sum, 64^3, against the point model alone at c = 1
(check_0d.json, the calibration's own trajectory):

| | So | An | Vei | Fn | Pg | sum | film c | thickness um |
|---|---|---|---|---|---|---|---|---|
| film, T* = 0 (day 1) | 0.613 | 0.006 | 0.065 | 0.053 | 0.002 | 0.739 | 0.46 | 16 |
| film, T* = 1 (day 21) | 0.386 | 0.048 | 0.109 | 0.013 | 0.007 | 0.562 | 0.10 | 45 |
| point model alone, day 21 | 0.141 | 0.212 | 0.154 | 0.002 | 0.002 | 0.511 | 1 | - |

Reading: in the film the nutrient is 0.1 (depletion length 10 um, film 45 um), so the calibrated community changes much
less than the point model alone at c = 1: S. oralis stays dominant (living share 0.69 against 0.28) and A. naeslundii,
which takes over in the calibration, stays at 0.09 (0.42 alone). The film's position in the field changes the predicted
composition. The two time links disagree: with the calibration's link (T* = 1 = 20 days) the film grows 29 um in 20 days,
slower than in situ plaque (35 um in 6 days, see above). 128^3 below when it finishes.

128^3 (22 min): film at T* = 1 phibar (So, An, Vei, Fn, Pg) = 0.398, 0.041, 0.091, 0.012, 0.006 (sum 0.548; living shares
0.73, 0.08, 0.17, 0.02, 0.01), c 0.13, thickness 49 um; 64^3: shares 0.69, 0.09, 0.19, 0.02, 0.01, c 0.10, 45 um. The finding
(So stays dominant, An does not take over) holds on both grids; the shares change by up to 4 points between them.

## fifa part 11: 192^3 (layered; mixed without antibiotic) and 256^3 (mixed, antibiotic), 8 Nov

Table and GCIs in KEIO_CHAIN_LOG.md (8 Nov). In short: the film stress converges (mixed, antibiotic, 2 h_f: 12.57, 12.67,
12.69 mPa on 128/192/256^3, GCI 0.1 %; layered 0.3 %; no antibiotic 0.3 %); the living fraction converges in the layered
film (0.307, 0.321, 0.339 on 64-192^3; extrapolated 0.346, GCI 2.3 %) but not in the mixed film with antibiotic: 0.201,
0.284, 0.320, 0.365 on 64-256^3, the increments growing per cell size (the 128/192/256^3 "GCI 14 %" is not meaningful).
Without antibiotic the mixed film's living fraction (0.81, 0.83, 0.84) does not reach the asymptotic range either.
Paper: stress converged, layered living fraction converged, mixed living fraction only bounded below (0.37).

## Stress with measured stiffness (13 Nov 2026)

The mechanics is one way and the strains are small, so at a given growth field the stress is proportional to E.
The model uses E = 10 Pa (Klempt 2024); AFM on oral biofilms gives 0.6 to 41 kPa (Pattem et al. 2018,
`references/EXPERIMENTAL_DATA.md`). The table rescales the computed von Mises stresses; nothing is recomputed, and the growth
field is assumed unchanged (no feedback of the stress on growth).

| example | E = 10 Pa (computed) | E = 0.6 kPa | E = 14 kPa | E = 41 kPa |
|---|---|---|---|---|
| cube, separate seed (extrapolated) | 0.6 mPa | 36 mPa | 0.85 Pa | 2.5 Pa |
| cube, mixed seed (superbee, extrapolated) | 0.38 mPa | 23 mPa | 0.53 Pa | 1.6 Pa |
| cervical window, biofilm mean | 13 mPa | 0.76 Pa | 18 Pa | 52 Pa |

So with measured stiffness the stresses are of the order of 0.02 to 50 Pa, still far below the contact pressures at which
confined microbes slow their growth (about 1 MPa for yeast, Delarue et al. 2016).
