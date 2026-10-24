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
