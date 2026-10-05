# Coupling status: what is done and what is not (5 Oct 2026)

Slide numbers refer to `slides_1005.tex` (5 Oct deck). Runs of 5 Oct:
`ansys_usermat/apdl/RUN_1005_IKMHIWI03.md` (record), numbers as JSON in
`ansys_usermat/apdl/results/2026-10-05_ansys/` (tables:
`summarize_runs_json.py`, figures: `figs_beta.py`).

One page for the question "what works for one and for two species". Details
are in the files named in each row.

The ANSYS model of every row below is the partner's 2 mm cube, 512 elements,
seed of 32 elements, nutrient held on the face y = -1 mm, only rigid-body
motion constrained (`assets/fig_model_setup.png`). Parameters: Klempt et al.
2024 Table 2 (k_alpha = 1e-3 per T*, E = 10 Pa, nu = 0.49) and, decided on
5 Oct, beta = 0.02 mm^2/T* for the diffusion of phi (Table 2's beta = 2
converted to the 2 mm cube, an assumption; the element's example value 1e-4
only as a sensitivity case). The two-species decks keep the element's
example stiffness (1000 MPa), so their stresses are not compared with the
one-species ones; the composition does not depend on the stiffness.

## One species

| | what | status | where |
|---|---|---|---|
| 1 | growth law F_g = alpha I -> stress, as an ANSYS user material | done: constrained cube equals the closed form, free growth gives no stress (~1e-10) | slide 4, `ansys_usermat/apdl/README.md` |
| 2 | Eq. 36 integrated in the material call of the partner's element | done: equals the exact solution sinh/cosh in ANSYS (alpha to 3e-11) | Check 1 (slide 10) |
| 3 | whole model, Eq. 36 against the element's own growth variable | done: same stress pattern (cosine 1.0000), ratio 0.455 explained (1/2 x 10/11) | Check 2 (slide 11) |
| 4 | seeded element and its neighbours | done: alpha - 1 = k_alpha phi t; neighbours carry about 25x its von Mises stress | Check 3 (slide 12) |
| 5 | sign of the neighbours' stress (the paper's ring of tension) | **done (5 Oct)**: at beta = 0.02 every element of the first layer around the seed is in tension, mean stress +6.0e-5 Pa, 8^3 and 16^3 within 0.5 % (dt 0.025). At beta = 1e-4 56 % of the first layer (8^3) | `RUN_1005_IKMHIWI03.md`, `assets/fig_beta_seed_stress.png` |
| 6 | stiffness line `bio1 + bio1` (typo, should be `bio1 + bio2`) | **done**: fixed in every executable since 1 Oct (the paper-value runs included); confirmed numerically on 5 Oct | `RUN_1005_IKMHIWI03.md` item 6 |
| 7 | biofilm spreading towards the nutrient | **not shown in ANSYS**: the element's front term is inactive (reading slip); with the slip fixed in a test copy the front is slow and its direction is undefined where the nutrient is uniform. Recorded as a limitation. At beta = 0.02 phi spreads around the seed by diffusion (not towards the nutrient) | `ansys_usermat/apdl/FRONT_TERM_FIX.md` |
| 8 | mesh and time step | **done (5 Oct)**: at beta = 0.02 the seed mean stress agrees between 8^3 and 16^3 within 0.4 % (dt 0.025), the seed von Mises does not (+22 %); dt 0.1/0.05/0.025 converge at first order (dt 0.025 within 1 % of the extrapolated mean stress). At beta = 1e-4 nothing converges: the diffusion length (~0.01 mm) is below every mesh and the alpha step at the seed surface grows with refinement. 24^3 at beta = 0.02 pending (5 Oct evening); 32^3 does not fit in memory | `RUN_1005_IKMHIWI03.md` |
| 9 | nu = 0.45 instead of 0.49 | done: at beta = 0.02 the seed stresses change by 2-3 % (the thesis' 27 % is the beta = 1e-4 value) | `RUN_1005_IKMHIWI03.md` |

So for one species the coupling itself is done and verified, and with
beta = 0.02 the seed's mean stress and the ring of tension are mesh-independent.
Growth happens where the biofilm already is (the seed) and where phi diffuses
to; it does not spread towards the nutrient in the ANSYS runs.

## Two species

| | what | status | where |
|---|---|---|---|
| 1 | point model of Klempt et al. 2026, two species | done: cases 3 and 6 reproduce the paper (case 3: 0.5945 / 0.3926 against 0.60 / 0.40) | `JAXFEM/klempt2026_reproduction.py` |
| 2 | ANSYS calls the point model at every Gauss point (Fortran -> Python) | done: every call replayed in Python, identical bit for bit | slide 9, `ansys_usermat/apdl/SPLIT_COUPLING.md` |
| 3 | coupled ANSYS run: amount from the field, composition from the point model | done: cases 3 and 6, now with the point model in every element (5 Oct: 256/256 seed elements on 16^3) | slide 15, `RUN_1005_IKMHIWI03.md` |
| 4 | the two assumptions s = 0.15 and phi_cap = 0.9 | done: phi_cap 0.85-0.95 moves the seed mean share by up to 8 %; **s decides the share at T* = 1** (beta = 0.02, consumption 6: seed mean 0.49 / 0.32 / 0.07 for s = 0.05 / 0.15 / 0.5) | slide 16, `RUN_1005_IKMHIWI03.md` |
| 5 | composition where the biofilm spreads | **partly (5 Oct)**: at beta = 0.02 the point model runs in 304 of 512 elements; behind the diffusing edge the share stays near its start value (0.2 / 0.5 / 0.7 -> 0.21 / 0.50 / 0.70 two layers out), and the same parameter set the seed's start. A separate start value for late points (prop(35)) is in the fragment and tested; its ANSYS runs are pending (5 Oct evening) | `RUN_1005_IKMHIWI03.md`, `assets/fig_beta_share_front.png` |
| 6 | stiffness with two species (`bio1 + bio1` typo) | **done**: fixed since 1 Oct | `RUN_1005_IKMHIWI03.md` item 6 |
| 7 | composition acting on the stress (two-way coupling) | not in the thesis by choice: Klempt et al. 2024 drop the stress feedback (Eq. 33); outlook for Keio | `KLEMPT2024_REPRODUCTION.md` sec. 12 |
| 8 | the field's nutrient in the point model (two-way step 1) | **done in ANSYS (5 Oct)**: the Gauss-point nutrient goes to the point model as c* = c*_0 min(c/c_ref, 1) (an assumption). Case 6: share 0.025-0.142 in the seed (beta = 1e-4, consumption 4, 8^3; 16^3: 0.024-0.155); beta = 0.02, consumption 6: 0.12-0.46 (8^3), 0.07-0.47 (16^3), low share near the nutrient face. Case 3 hardly changes (0.62-0.66). Consumption values are example inputs; the element solves Eq. 35 without c-dot (zero order), so c < 0 above consumption ~6 | `RUN_1005_IKMHIWI03.md`, `ROADMAP_TWO_WAY.md` |

So for two species the coupling is complete and runs correctly in ANSYS, with
a composition that varies in space through the local nutrient. The share at
T* = 1 depends strongly on the assumed clock s; behind the diffusing edge it
depends on the assumed start value.

## Open decisions (5 Oct)

- Which s to report as the main case (the share at T* = 1 follows s·T*).
- How to treat the start value behind the edge (prop(35) runs pending; Felix
  asked by email on 5 Oct how new points start in Klempt et al. 2026).
- beta: conversion of Table 2's beta = 2 to the 2 mm cube (Felix and Oliver
  asked on 5 Oct).

## More than two species

Not in the thesis (continued at Keio). The parts exist: the five-species
point model (TMCMC-calibrated, c* = 25), the same bridge, and five-species
runs on this repository's own cylinder test decks. With the scheme used for
two species (amount from the field, composition from the point model) the
field stays a single phi, so only the point-model state and its parameters
would need widening.
