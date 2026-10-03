# Coupling status: what is done and what is not (3 Oct 2026)

Slide numbers refer to `slides_1005.tex` (5 Oct deck).

One page for the question "what works for one and for two species". Details
are in the files named in each row.

The ANSYS model of every row below is the partner's 2 mm cube, 512 elements,
seed of 32 elements, nutrient held on the face y = -1 mm, only rigid-body
motion constrained (`assets/fig_model_setup.png`). Parameters: Klempt et al.
2024 Table 2 (k_alpha = 1e-3 per T*, E = 10 Pa, nu = 0.49).

## One species

| | what | status | where |
|---|---|---|---|
| 1 | growth law F_g = alpha I -> stress, as an ANSYS user material | done: constrained cube equals the closed form, free growth gives no stress (~1e-10) | slide 4, `ansys_usermat/apdl/README.md` |
| 2 | Eq. 36 integrated in the material call of the partner's element | done: equals the exact solution sinh/cosh in ANSYS (alpha to 3e-11) | Check 1 (slide 10) |
| 3 | whole model, Eq. 36 against the element's own growth variable | done: same stress pattern (cosine 1.0000), ratio 0.455 explained (1/2 x 10/11) | Check 2 (slide 11) |
| 4 | seeded element and its neighbours | done: alpha - 1 = k_alpha phi t; neighbours carry about 25x its von Mises stress | Check 3 (slide 12) |
| 5 | sign of the neighbours' stress (the paper's ring of tension) | **open**: read on IKMHIWI03 on Monday | `PAPER_CHECK_KLEMPT2024.md` |
| 6 | stiffness line `bio1 + bio1` (typo, should be `bio1 + bio2`) | **open**: fix on Monday, then check whether the one-species stresses change | `PAPER_CHECK_KLEMPT2024.md` |
| 7 | biofilm spreading towards the nutrient | **not shown in ANSYS**: the element's front term is inactive (reading slip); with the slip fixed in a test copy the front is slow and its direction is undefined where the nutrient is uniform. Recorded as a limitation | `ansys_usermat/apdl/FRONT_TERM_FIX.md` |

So for one species the coupling itself is done and verified. Growth happens
where the biofilm already is (the seed); it does not spread in the ANSYS runs.

## Two species

| | what | status | where |
|---|---|---|---|
| 1 | point model of Klempt et al. 2026, two species | done: cases 3 and 6 reproduce the paper (case 3: 0.5945 / 0.3926 against 0.60 / 0.40) | `JAXFEM/klempt2026_reproduction.py` |
| 2 | ANSYS calls the point model at every Gauss point (Fortran -> Python) | done: every call replayed in Python, identical bit for bit | slide 9, `ansys_usermat/apdl/SPLIT_COUPLING.md` |
| 3 | coupled ANSYS run: amount from the field, composition from the point model | done in the seed region: cases 3 and 6, 0 errors, all trace checks pass (1 Oct) | slide 15 |
| 4 | the two assumptions s = 0.15 and phi_cap = 0.9 | done: sensitivity studies (case 3: 0.664-0.668 for phi_cap 0.85-0.95) | slide 16, `ansys_usermat/composition_figs.py` |
| 5 | composition where the biofilm grows (the front) | **not shown in ANSYS**, for the same reason as row 7 above. Shown in Python on the paper's 20 um cube: starting new points from the seed or from their neighbour changes the mean composition by at most 0.001 | slide 17, `ansys_usermat/composition_transport_check.py` |
| 6 | stiffness with two species (`bio1 + bio1` typo) | **open**: Monday | `PAPER_CHECK_KLEMPT2024.md` |
| 7 | composition acting on the stress (two-way coupling) | not in the thesis by choice: Klempt et al. 2024 drop the stress feedback (Eq. 33); outlook for Keio | `KLEMPT2024_REPRODUCTION.md` sec. 12 |

So for two species the coupling is complete and runs correctly in ANSYS. The
results that can be shown in ANSYS are those in the packed seed; the
composition along a spreading front cannot be shown there until the
element's front term works.

## More than two species

Not in the thesis (continued at Keio). The parts exist: the five-species
point model (TMCMC-calibrated, c* = 25), the same bridge, and five-species
runs on this repository's own cylinder test decks. With the scheme used for
two species (amount from the field, composition from the point model) the
field stays a single phi, so only the point-model state and its parameters
would need widening.
