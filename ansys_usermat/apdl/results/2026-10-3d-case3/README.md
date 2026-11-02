# Coupled runs with the inputs of Klempt et al. 2026, Table 2, case 3 (kept apart from the main results)

The user (21 Oct 2026): "align everything with Klempt 2026 Table 2 case 3; try both 2 and 4 species; save as a
separate result". All runs: paper values of `references/PARAMETER_TABLES.md` (r, d, g of Klempt 2024 Table 2),
quasi-static fields, Python `coupled_living_consumption.py`.

| tag | meaning |
|---|---|
| `_csine` | held nutrient 0.5 + 0.5 sin(500 t), t = s T* the point-model time, i.e. c* = 50 + 50 sin(500 t) (case 3) |
| `_tab0` | antibiotic from the start, alpha*_0 = 10: alpha* = 10 a (case 3: alpha* = 10 constant) |
| `_k26c3` | the four-species point model of case 3: A = 1/2 [[1,5,5,5],[5,1,3,3],[5,3,1,2],[5,3,2,1]], eta = (0.8, 1, 1.5, 2), b = (0.4, 0.3, 0.2, 0.1), start phi 0.02 each, 1500 steps of 1e-4 (s = 0.15); species 1, 2 = "A" (seed), 3, 4 = "B" (corners); "share of A" = (phi_1 + phi_2) / sum phi_i |
| without `_k26c3` | two species of Fritsch case I with the case-3 inputs (s = 0.1) |
| `_dtc0.005`, `_dtc0.004` | outer step: the sine has a period of 0.126 T* (2 species) and 0.084 T* (4 species) |

Check: with c = a = 1 everywhere the seed point model equals the point model alone to 1.5e-15 (4 species, outer
step 0.004) and 1.7e-14 (2 species); the reference uses a converged Newton (6 fixed iterations are not converged for
four species: 5.6e-6). Assumption: t in sin(500 t) is the point-model time in units of Time (dt = 1e-4), as the
paper's text states the step size in Time; the paper's "t > 500" of case 4 counts steps.

**24 Oct 2026 (A-77):** the four-species matrix A is used without the printed 1/2 of Eq. 19 (only then are the
paper's Figs. 10 and 11 reproduced, 0D check `checks_a38_a39/k26_point_check_a77.py`). Files with `_k26c3_` (no `A1`)
used the 1/2 and are superseded; current: `_k26c3A1_` (case 3) and `_k26c4_` (case 4, the antibiotic case).

## Two species to four species (case 4, `_k26c4`): what changes (2 Nov 2026)

Only the point model and what follows from it; the field equations (Klempt 2024 Eqs. 34-36), the constants, the grid,
the boundary conditions, the antibiotic from T* = 1/3 (alpha*_0 = 100) and c*_0 = 100 are those of the two-species runs.

- Point-model case: Fritsch et al. 2025 case I (two species) -> Klempt et al. 2026 Table 2 case 4 (A as above without
  the printed 1/2, eta = (0.8, 1, 1.5, 2), b = (10, 2, 1, 0.01), start phi 0.02 each).
- Time link s = 0.1 -> 0.15 (T* = 1 spans the 1500 steps of the case), so the inner point-model steps per outer step
  0.0025 grow from 25 to 38.
- Start: separate = species 1, 2 in the seed ("A"), species 3, 4 in the edge blocks ("B"); mixed (`--mix-equal`) = all
  four in equal amounts (0.02 each) in the seed.
- Unchanged in form: uptake g phi_3D (sum_i phibar_i) c (the sum now over four species), the front term
  -r c/(k+c) n_c . grad phi_3D, the growth alpha-dot = k_alpha phi_3D, the stress solve. 12 state values per point in
  both cases.
- Extra output only: shares and phibar_i per species (`sp{i}` snapshots, `share_sp_*`, `phibar_sp_*` histories).
- Cost at 64^3 (quarter domain, outer step 0.0025, antibiotic): about 350-470 s against 170-350 s for two species
  (1.5-2x, from the 38 instead of 25 inner steps); memory about the same.
- End values (T* = 1): seed living 0.76 / 0.80 (separate / mixed), phibar_4 over the biofilm 0.75 / 0.79 against 0.93
  for the point model alone (c ~ 0.6 in the seed), seed von Mises 6.5e-4 / 3.5e-4 Pa. Deck: `meetings/2026-11_coupled_examples/deck_k26c4_en.js`.
