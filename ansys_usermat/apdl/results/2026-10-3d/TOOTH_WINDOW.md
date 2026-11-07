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
