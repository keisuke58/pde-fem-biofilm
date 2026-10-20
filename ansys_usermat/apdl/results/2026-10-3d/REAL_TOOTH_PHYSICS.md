# Nutrient field and homeostatic-pressure gate on the real molar (3D demonstration, 20 Oct 2026)

Script `ansys_usermat/apdl/real_tooth_model.py` (the seed model of `anim_real_tooth_3d.py`
with a nutrient field and the pressure gate added, each on or off); numbers
`real_tooth_physics.json` (this directory, one entry per run with the z profiles);
videos `assets/anim_real_tooth_gate_3d.mp4` (panels phi, alpha - 1, von Mises) and
`assets/anim_real_tooth_nutrient_gate_3d.mp4` (c, alpha - 1, von Mises), 2400 x 912, 16 fps,
90 frames, with their last frames as `.png`; the last frame of the nutrient-only run
`assets/anim_real_tooth_nutrient_3d.png`; figure `assets/fig_real_tooth_physics.pdf` / `.png`;
caches `figdata/real_tooth_<tag>.npz` (git-ignored, rebuilt by the script when missing; the
base run reads `figdata/anim_real_tooth_3d.npz`). This is a demonstration on the molar
FDI 46 of `REAL_TOOTH.md` (real size in a 16 mm cube, CC BY 4.0 geometry), not a thesis
result.

## Equations as implemented

Field on 96^3 cells (h = 0.167 mm), classical RK4, dt = 0.0167 T* (one step per frame, 66
frames to T* = 1.1), zero flux at the box walls and at the tooth surface:

- biofilm state variable (Eq. 34 of Klempt et al. 2024 without the front term, as in all the
  cube demonstrations):
  `d phi/dt = beta lap(phi) + k_alpha alpha - penalty (max(phi - 1, 0) + min(phi, 0))`;
- nutrient (Eq. 35, quasi-static as in `NUTRIENT_3D.md`; consumption `g phi c`, which
  Dr. Soleimani confirmed on 6 Oct 2026 as the authors' implementation of Eq. 24):
  `0 = d lap(c) - g phi c` in the open cells, c = 1 on the five outer faces of the box except
  the clamped bottom (the saliva around the crown), zero flux at the bottom face and at the
  tooth. Only g/d enters. Solved once per RK4 step (sparse 7-point matrix, pyamg smoothed
  aggregation with CG, relative residual below 1e-9, warm-started: 7 to 11 iterations, about
  6 s per solve);
- local expansion (Eq. 36 with the factor c and the gate g):
  `d alpha/dt = k_alpha phi c g`, alpha(0) = 1, F_g = alpha I. Eq. 36 as printed is
  `d alpha/dt = k_alpha phi`; the factor c is an assumption of these demonstrations (a linear
  nutrient dependence of the local growth), as in `NUTRIENT_3D.md`. c = 1 or g = 1 is held
  where the nutrient or the gate is off;
- gate (the law case of `anim_scenarios_3d.py`, the homeostatic pressure law of P3):
  `g = min(1, max(0, 1 - p/p_h))`, p = -tr(sigma)/3 (compression positive), from the stress at
  the start of the frame and held over the frame, per stress element (the 0.5 mm element
  value is applied to its 3^3 field cells). With the gate the stress is solved every frame
  (67 solves) instead of every third frame. The 32^3 stress mesh was fast enough for that
  (12 s per solve with one BLAS thread), so the stress mesh is the same in all runs;
- mechanics as in `REAL_TOOTH.md`: small strain, hex8 B-bar on 32^3 (h = 0.5 mm, the fields
  block-averaged), E = 10 Pa (phi^2 + 1e-3), tooth 1e4 Pa without growth, nu = 0.49,
  eigenstrain (alpha - 1) I, bottom face z = -8 mm clamped, pyamg + CG, tolerance 1e-8.

## Parameters

| quantity | value | source |
|---|---|---|
| beta | 0.08 mm^2/T* | Table 2 converted to the 2 mm cube (0.02 mm^2/T*) and scaled x 4 to the 0.5 mm film: an assumption, explained in `REAL_TOOTH.md` |
| k_alpha | 1e-3 /T* | Table 2 |
| penalty | 5 | the partner's deck, as in the other demonstrations |
| g/d, main runs | 10 /mm^2 | ten times the Table 2 ratio (g = 1e8 /T*, d = 1e10 um^2/T*, converted like beta to 1 /mm^2, `NUTRIENT_3D.md`); chosen by me so that the penetration depth (d/g)^1/2 = 0.32 mm is below the film thickness and a gradient shows across the 0.5 mm film |
| g/d, check run | 1 /mm^2 | the Table 2 ratio (penetration depth 1 mm); numbers only |
| p_0 | 3.705 mPa | mean pressure in the film (phi > 0.5 on the stress mesh, 37 elements) at the end of the base run |
| p_h | 1.853 mPa = 0.5 p_0 | assumption. The threshold study of `p3/uel/README.md` ("Threshold p_h = p_0") placed the crossing of the band-to-core stress ratio near p_h = p_0 with p_0 the band mean; here I take the film mean and half of it, so that the gate acts within the run |
| E, nu | 10 Pa, 0.49; tooth 1e4 Pa | Table 2; the tooth stiffness is an assumption (rigid against the film) |
| stiffness floor | 1e-3 | assumption, as in the other demonstrations |
| T* | 1.1 | as in the other demonstrations |

Film and band for the numbers: film = cells with phi > 0.5 on the field grid (alpha - 1, c;
1857 cells) and elements with phi > 0.5 on the stress mesh (stress, pressure; 37 elements,
the film is one element thick); the wider film phi > 0.3 (130 elements) is used for the von
Mises profile of the figure and given in the json. Cervical band = film with z in
[-2.1, -1.1] mm (the cemento-enamel junction is near z = -1.6 mm; 18 elements in two element
layers), rest = the film outside it (19 elements). I report means and percentiles, not the
pointwise maximum: `WALL_CONVERGENCE.md` showed that the maximum at clamped or free edges
grows with the stress mesh while the band means converge. With 37 and 18 elements the
percentiles are coarse; the means are the numbers to read.

## Results at T* = 1.1

Stresses in mPa (1e-3 Pa), alpha - 1 in 1e-3. The film volume (sum of phi) is 17.24 mm^3 in
every run (the phi equation does not depend on c or g); 1857 film cells, 130 / 37 elements.

| run | c in the film, min / mean / max | film mean alpha - 1 | film mean vM | film 95th pct vM | film mean p | p / p_h | share of the film with p >= p_h | gate factor g, film mean |
|---|---|---|---|---|---|---|---|---|
| no nutrient, no gate (base) | 1 | 0.858 | 3.81 | 5.87 | 3.71 | 2.00 | 0.57 | 1 |
| nutrient, g/d = 10 /mm^2 | 0.067 / 0.135 / 0.237 | 0.116 | 0.54 | 0.77 | 0.40 | 0.21 | 0.00 | 1 |
| gate | 1 | 0.521 | 1.68 | 2.16 | 1.18 | 0.64 | 0.00 | 0.36 |
| nutrient + gate | 0.067 / 0.135 / 0.237 | 0.104 | 0.46 | 0.60 | 0.33 | 0.18 | 0.00 | 0.82 |
| nutrient, Table 2 ratio 1 /mm^2 (numbers only) | 0.624 / 0.692 / 0.782 | 0.580 | 2.58 | 3.89 | 2.42 | 1.31 | 0.35 | 1 |

Cervical band against the rest of the film (von Mises, mPa):

| run | band mean | band 95th pct | rest mean | rest 95th pct | ratio band / rest, mean | ratio, 95th pct | band mean p | rest mean p | band mean alpha - 1 | rest mean alpha - 1 |
|---|---|---|---|---|---|---|---|---|---|---|
| base | 3.57 | 5.86 | 4.03 | 5.86 | 0.89 | 1.00 | 2.90 | 4.47 | 0.859 | 0.855 |
| nutrient | 0.54 | 0.72 | 0.55 | 0.79 | 0.98 | 0.91 | 0.31 | 0.48 | 0.115 | 0.118 |
| gate | 1.70 | 2.17 | 1.66 | 2.06 | 1.02 | 1.05 | 1.12 | 1.24 | 0.536 | 0.500 |
| nutrient + gate | 0.46 | 0.57 | 0.46 | 0.62 | 1.00 | 0.91 | 0.28 | 0.38 | 0.104 | 0.104 |
| nutrient, Table 2 ratio | 2.44 | 3.79 | 2.71 | 3.95 | 0.90 | 0.96 | 1.90 | 2.91 | 0.581 | 0.580 |

Seed means (the initial seed cells, 2954): alpha - 1 = 0.757 (base), 0.125 (nutrient), 0.489
(gate), 0.113 (nutrient + gate), 0.526 (Table 2 ratio); c = 0.170 in the seed at g/d = 10 and
0.713 at the Table 2 ratio.

Von Mises by height (film phi > 0.3, element layers z = -2.75, -2.25, -1.75, -1.25 mm,
panel (a) of the figure): base 1.68, 2.65, 2.64, 1.43; nutrient 0.29, 0.38, 0.40, 0.27;
gate 0.85, 1.20, 1.30, 0.85; nutrient + gate 0.25, 0.33, 0.35, 0.24; Table 2 ratio 1.17,
1.80, 1.80, 1.00 mPa. The profile keeps its shape (largest in the middle two layers, lower
at the top and bottom edges of the band) in every run; the nutrient and the gate scale it.
c by height (panel (c)): 0.12 to 0.20 at g/d = 10 (lowest in the middle of the band, where
the film is thickest and farthest from the open faces), 0.68 to 0.74 at the Table 2 ratio.

Wall times (one BLAS thread, two runs at a time on the shared 4 cores): gate 19 s field +
701 s stress (67 solves); nutrient 21 s + 390 s nutrient (67 solves) + 254 s stress (23
solves); nutrient + gate 21 + 393 + 701 s; Table 2 ratio 21 + 373 + 12 s (stress at the end
only); videos 73 and 103 s.

## What the nutrient and the gate change on the real tooth

The nutrient supplied by the saliva around the crown is almost used up before it reaches the
film interior at g/d = 10 /mm^2: c falls from 1 at the box faces to 0.24 at the outer surface
of the film and 0.07 at the tooth, 0.135 on average, so the growth and with it the stress
drop to about one seventh of the base run (alpha - 1 from 0.86 to 0.12 x 1e-3, film mean von
Mises from 3.8 to 0.5 mPa) and the film pressure (0.4 mPa) no longer reaches the threshold
p_h. At the Table 2 ratio the nutrient reaches the whole film at c = 0.62 to 0.78 and growth
and stress drop by about a third (0.68 x), nearly uniformly, so the stress pattern of the
base run stays. The gate alone throttles the growth once the film pressure rises towards p_h:
the growth factor g falls to 0.36 on average at the end, the film mean alpha - 1 to 0.61 x the
base value and the film mean von Mises to 0.44 x (1.7 mPa), while the pressure levels off at
0.64 p_h (the law is linear in p, so p approaches p_h from below and the share of the film
with p >= p_h is 0 by construction, against 0.57 in the base run). The gate acts most where
the pressure is highest, which is the rest of the film rather than the cervical band (base
run: 4.5 against 2.9 mPa), so it evens the stress out: the band-to-rest ratio of the mean
von Mises goes from 0.89 to 1.02, and the band keeps slightly more growth than the rest
(0.54 against 0.50 x 1e-3). With both, the nutrient dominates: the pressure stays at 0.18 p_h,
the gate hardly acts (g = 0.82 on average) and the result is within 12 % of the nutrient run
alone (alpha - 1 0.104 against 0.116 x 1e-3, von Mises 0.46 against 0.54 mPa).

## Limitations

- The film on the 0.5 mm stress mesh is one element thick, 37 elements at phi > 0.5 and 18
  in the cervical band, so the band statistics rest on few elements and the 95th percentiles
  are coarse. The von Mises pattern on the film surface in the videos shows the element grid,
  as in `REAL_TOOTH.md`.
- The gate is applied per stress element and held over a frame, so the growth of the gate
  runs carries the element pattern: the alpha - 1 panel of the gate video shows vertical
  stripes of the 32^3 mesh. A finer stress mesh or a smoothed gate would remove them; the
  means are not affected.
- The gate is explicit (pressure from the start of the frame): p overshoots by at most one
  frame of growth, which is small here because p stays below p_h. T* = 1.1 is not long enough
  for p to settle at p_h; the gate run is still approaching it (0.64 p_h at the end).
- p_h = 0.5 p_0 with p_0 the film mean of the base run is a choice made so that the gate
  acts within the run; the threshold study of `p3/uel/README.md` used the cervical band mean.
  A different p_h scales the gate's effect.
- The factor c in d alpha/dt is an assumption (Eq. 36 has no c; in the paper the nutrient
  acts through the front term of Eq. 34, which these demonstrations do not carry). g/d = 10
  is ten times the Table 2 ratio and was chosen to show a gradient across the film; at the
  Table 2 ratio the limitation is a nearly uniform factor of 0.68.
- c is quasi-static and consumed only by phi (no saliva flow, no clearance, no supply from
  the gum side); c = 1 on the box faces places the saliva 1.7 mm or more from the film, and
  the tooth is impermeable. The share of the drop of c that happens in the gap between the
  box faces and the film is part of the result (as in `NUTRIENT_3D.md`).
- beta is scaled x 4 (see `REAL_TOOTH.md`); the tooth is rigid and not grown; nothing below
  the gum line is modelled; the space outside the film has the floor stiffness, not a fluid.
- The film volume does not depend on c or g because the phi equation of these demonstrations
  has no front term; only the local expansion alpha and the stress respond.

## How to rerun

```
python ansys_usermat/apdl/real_tooth_model.py gate                 # run + video, about 13 min
python ansys_usermat/apdl/real_tooth_model.py nutrient_gate        # about 20 min
python ansys_usermat/apdl/real_tooth_model.py nutrient --png-only  # about 11 min
python ansys_usermat/apdl/real_tooth_model.py nutrient_t2 --no-video
python ansys_usermat/apdl/real_tooth_model.py base --no-video      # numbers from the base cache
python ansys_usermat/apdl/real_tooth_model.py --summary            # json + figure from the caches
```
Each run skips itself when its cache exists and merges its numbers into
`real_tooth_physics.json` as soon as it finishes; `--summary` rebuilds the json and the
figure from all caches. Needs pyvista (off screen), imageio-ffmpeg, scipy and pyamg.
