# Seed demonstrations with a nutrient field (3D, 2 mm cube)

Demonstrations with the Chapter 4 seed model, not thesis results. Script:
`ansys_usermat/apdl/anim_nutrient_3d.py`; numbers: `nutrient_<case>.json` in this
directory; videos and last frames: `assets/anim_nutrient_<case>_3d.mp4/.png`;
cached fields: `figdata/anim_nutrient_<case>_3d.npz`. Date: 20 Oct 2026.

I added the nutrient concentration c of Klempt et al. (2024) to the cube
demonstrations of `anim_scenarios_3d.py`, so the growth is no longer uniform in
the biofilm but limited by the nutrient that diffuses in from the top face and
is consumed by the biofilm.

## Equations as implemented

Field on 48^3 cells of the cube [-1, 1]^3 mm, classical RK4 in time, T* in
[0, 1.1]:

- biofilm state variable (Eq. 34 of Klempt et al. 2024 without the front term,
  as in all the cube demonstrations):
  `d phi/dt = beta lap(phi) + k_alpha alpha - penalty (max(phi - 1, 0) + min(phi, 0))`,
  zero flux on all six faces;
- nutrient (Eq. 35, quasi-static, i.e. the paper's time homogenisation with
  eta_c = 1e-10; consumption `g phi c`, which Dr. Soleimani confirmed on 6 Oct
  2026 as the authors' implementation of Eq. 24):
  `0 = d lap(c) - g phi c`, c = 1 on the top face z = +1 (saliva side), zero
  flux on the other five faces. Only the ratio g/d enters;
- local expansion (Eq. 36 with the factor c):
  `d alpha/dt = k_alpha phi c`, alpha(0) = 1, F_g = alpha I.
  Eq. 36 as printed is `d alpha/dt = k_alpha phi`; in the paper the nutrient
  acts on the growth through the front term `r c/(k + c)` of Eq. 34, which
  these demonstrations do not carry. The factor c is therefore an assumption of
  this demonstration (a linear nutrient dependence of the local growth). The
  run with c = 1 held everywhere is the model of `anim_scenarios_3d.py`.
- mechanics as before: small strain, hex8 B-bar on 16^3 (fields block-averaged),
  E = 10 Pa (phi^2 + 1e-3), nu = 0.49, eigenstrain (alpha - 1) I, bottom face
  z = -1 clamped, solved once per frame (88 frames).

Each case is run twice, with the nutrient and with c = 1 held (no limitation),
so that the effect can be given in numbers.

## Parameters

| quantity | value | source |
|---|---|---|
| beta | 0.02 mm^2/T* | Table 2 (2 um^2/T*) converted to the 2 mm cube, KLEMPT2024_REPRODUCTION.md section 8 |
| k_alpha | 1e-3 /T* | Table 2 |
| penalty | 5 | the partner's deck (PENALTY1), as in the other demonstrations |
| g/d | 1 /mm^2 | Table 2: g = 1e8 /T*, d = 1e10 um^2/T*; converted like beta (lengths x 100, so d = 1e8 mm^2/T*), g L^2/d is preserved (`klempt2024_cases.deck_values`) |
| g/d, cases `_g10` | 10 /mm^2 | g x 10 (1e9 /T*), between the paper's two Fig. 7 values 1e8 and 1e10; chosen by me to show a gradient across the layer thickness |
| E, nu | 10 Pa, 0.49 | Table 2 |
| stiffness floor | 1e-3 | assumption, as in the other demonstrations |
| T* | 1.1 | as in the other demonstrations |

The penetration depth (d/g)^1/2 is 1 mm at the Table 2 ratio and 0.32 mm at
g x 10. With the supply on the top face, most of the drop of c happens in the
void between the top face and the biofilm (1.5 to 1.75 mm of pure diffusion
with the flux the biofilm consumes), so at the Table 2 ratio c is nearly
uniform inside a 0.25 mm layer and the gradient across the layer thickness
only shows at g x 10. Initial condition: c is quasi-static, so the c field at
T* = 0 is the steady solution for the initial seed (no initial value is
needed).

Cases:

- `corner`: seed x < 0, y < 0, z < -0.5 (1 x 1 x 0.5 mm) on the clamped bottom
  face, the corner seed of the `law` case; g/d = 1 /mm^2;
- `layer`: flat layer z < -0.75 (0.25 mm) over the whole bottom face; g/d = 1 /mm^2;
- `corner_g10`, `layer_g10`: the same seeds with g/d = 10 /mm^2.

## Stability check

- phi, alpha: RK4 with dt = 0.0125 T* (one step per frame); diffusion number
  beta dt/h^2 = 0.144 with h = 1/24 mm, below the RK4 limit 2.785/12 = 0.232
  for the 7-point stencil; penalty dt = 0.0625.
- c: an explicit step would need d dt/h^2 < 0.232, i.e. dt < 4.0e-12 T* with
  the converted d = 1e8 mm^2/T* (diffusion number 7e8 at the step used). This
  is why the paper treats c as quasi-static, and I do the same: the steady
  equation is solved once per RK4 step by conjugate gradients (Jacobi
  preconditioner, relative residual below 1e-10, at most 361 iterations cold
  and 48 warm-started for the layer). The reaction number q h^2 is 0.0017
  (0.017 at g x 10), so the consumption is resolved on the grid. The quasi-static
  solve has no step limit.
- the c solve is matrix-free; the threaded BLAS must be switched off for it
  (the script sets OMP/OPENBLAS/MKL threads to 1; with threads the 110k-vector
  dot products are 100 times slower).

Runtime per case: about 4 minutes for both runs (48^3 field, 16^3 stress, 88
frames) plus about 3 minutes for the video, on one core.

## Results at T* = 1.1 (seed = the initial seed cells of the stress mesh)

Values with the nutrient / with c = 1 held; stress in mPa (1e-3 Pa).

| case | c in the biofilm, min to max (mean) | seed mean alpha - 1 (1e-3) | seed mean vM | seed max vM | seed mean p | ratio alpha - 1, vM mean, vM max |
|---|---|---|---|---|---|---|
| corner, g/d = 1 | 0.71 to 0.81 (0.76) | 0.667 / 0.877 | 3.27 / 4.43 | 11.9 / 16.5 | 1.60 / 2.15 | 0.76, 0.74, 0.72 |
| layer, g/d = 1 | 0.685 to 0.698 (0.69) | 0.590 / 0.859 | 4.32 / 6.30 | 6.89 / 10.1 | 2.55 / 3.72 | 0.69, 0.69, 0.68 |
| corner, g/d = 10 | 0.12 to 0.34 (0.22) | 0.203 / 0.877 | 0.82 / 4.43 | 2.70 / 16.5 | 0.42 / 2.15 | 0.23, 0.18, 0.16 |
| layer, g/d = 10 | 0.155 to 0.186 (0.167) | 0.147 / 0.859 | 1.06 / 6.30 | 1.62 / 10.1 | 0.62 / 3.72 | 0.17, 0.17, 0.16 |

Profiles over the stress-mesh layers of the seed (z = -0.94, -0.81, -0.69,
-0.56 mm for the corner seed; -0.94, -0.81 mm for the layer), alpha - 1 in
1e-3, with the nutrient / with c = 1:

- corner, g/d = 1: 0.73, 0.72, 0.68, 0.53 / 0.98, 0.96, 0.89, 0.68 (c: 0.76, 0.76, 0.77, 0.79);
- corner, g/d = 10: 0.19, 0.20, 0.22, 0.20 / 0.98, 0.96, 0.89, 0.68 (c: 0.21, 0.22, 0.25, 0.29);
- layer, g/d = 1: 0.66, 0.52 / 0.96, 0.75 (c: 0.686, 0.694);
- layer, g/d = 10: 0.156, 0.138 / 0.96, 0.75 (c: 0.158, 0.176).

Location of the largest von Mises stress in the biofilm (mm): corner seed
(-0.81, -0.81, -0.94) with c = 1 and at g/d = 1, (-0.44, -0.44, -0.94) at
g/d = 10; layer (-0.06, -0.06, -0.94) in all runs, i.e. at the clamped bottom
face in every case.

## What changes

At the Table 2 ratio the nutrient reaches the whole seed at c = 0.7 to 0.8, so
the growth and the stress drop by 25 to 30 % almost uniformly and the stress
concentrates where it did before, at the clamped bottom face (the outer corner
of the corner seed, the centre of the layer). At ten times the consumption the
nutrient inside the seed falls to 0.12 to 0.34 and now decreases from the
upper and outer faces of the seed towards the clamped corner, which flattens
the growth over the seed height (alpha - 1 is 0.19 to 0.22 x 1e-3 in all four
layers, where without the nutrient it falls from 0.98 to 0.68 x 1e-3 from the
bottom to the top, because phi thins at the top by diffusion) and moves the
largest von Mises stress from the outer corner of the seed to the middle of
its footprint (-0.44, -0.44, -0.94 mm), nearer the faces that look into the
open void and so get more nutrient. The stress itself falls to about one sixth; the layer behaves the
same way, with c 11 % lower at its bottom than at its top at g x 10 and the
stress still largest at the clamped face.
