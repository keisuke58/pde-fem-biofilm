# Python runs of the agreed coupling (`coupled_living_consumption*.json`): which results to use

Updated 22 Oct 2026; the confirmed facts behind it: `../../CONFIRMED_2026-10-22.ja.md`.

**Update 25 Oct 2026 (thesis numbers):** new points inherit the neighbours' living composition (`_inhliv`, replaces
`_startpos`, A-82). Thesis Ch. 4 uses `_inhliv` at 64^3 with the outer step 0.025 for the example table and the grid
series `_dtc0.0025` (`*_front1_dtc0.0025_*inhliv*`, no prefix = 32^3, `nf64_`, `nf128_`; `nf192_` mixed_ab running)
for growth and stress (A-89; `../../CLOUD_TO_IKMHIWI03.ja.md` §2). The 22 Oct main line below is kept for history.

**Current main results: `_front1_phimin0.001_v2_fortran_ansysbc_startpos[_mixeq]_a100_gab1_gc1_tab0.333333_tvd_wgmirror_fast`**
at 128^3 (64^3 as the mesh check): paper values (front r of Klempt 2024 Table 2, d and g of Table 2, g l^2/d = 1,
alpha*_0 = 100 of Klempt 2026 case 4, antibiotic from T* = 1/3 with the nutrient's diffusion and uptake, phi_min = 0.001),
second-order upwind front term, wall gradient by the mirror (zero-flux) rule.

**Superseded on 22 Oct:** `_tvd_fast` and `_tvd` without `_wg`: the gradient of c (direction of the front) was one-sided
in the wall cells (stress changes 10 % at 16^3, 7 % at 32^3). `_artdiff` without `_wg`: the central front term used a
one-sided phi gradient in the wall cells, so wall blocks stayed at phi = 1.1 (A-38): **wrong, do not use**.

| tag | meaning |
|---|---|
| `_tvd` | second-order upwind (MUSCL, minmod) front term: **main** |
| `_artdiff` | as the ANSYS element: central front term + beta + r h/4; **only for comparisons with ANSYS** (at the ANSYS mesh, `--dtc` = the ANSYS increment, e.g. `_dtc0.005` at 8^3, `_dtc0.00125` at 16^3/24^3) |
| `_wgmirror`, `_wgzero` | normal component of the gradient in the wall cells: mirror (zero-flux ghost, main) or zero (the `--artdiff` default; the closer one to ANSYS at 8^3). No `_wg` = one-sided (wrong) |
| `_arth0.25` | `--artdiff-h`: r h/4 with the ANSYS element size at a finer Python grid (cloud test, 22 Oct) |
| `_fast` | quarter domain with mirror symmetry and AMG-CG (same results) |
| `_fortran` | the point model of the ANSYS runs through ctypes; without it the JAX version (agrees to 1e-16; IKMHIWI03 runs have no `_fortran`) |
| `_transc1e+08` / no `_transc` | dc/dt with d = 1e8 mm^2/T* / quasi-static: agree to < 1e-4 |
| `_transc1`, `_transc0.3`, `_transc3` | d = 1, 0.3, 3 mm^2/T* (Oliver's example deck value and around it): **slow-diffusion sensitivity only** |
| `_ansysbc` | boundary conditions of the ANSYS decks (layers held): use for comparisons with ANSYS |
| `_startpos` | new points start by position (the ANSYS rule); without it: from the neighbours' state |
| `_mixeq` | mixed start with equal amounts (the agreed example); without it in a mixed run: Fritsch's start (0.25, 0.30) |
| `_nf8` ... `_nf128` | mesh (default 32^3) |
| `_v2` | the v2 set-up (B in corner blocks, antibiotic on all outer faces) |
| `_a10`, no `_tab`, `_phimin0.1`, `_front0.1`, `_gc` other than 1 | earlier steps (alpha*_0 = 10, antibiotic from 0.5, phi_min 0.1, r/10, g l^2/d = 4): superseded |
| `_homeo<P>[f1|f2]` | Keio continuation, not in the thesis (27 Oct): the homeostatic-pressure growth law `--homeostatic P` (P = p_h/E_bio; Eq. 36 times g = min(1, max(0, 1 - p/p_h)), p = -tr(sigma)/3 of the cell from the stress solve of the start of the outer step, stress mesh 16^3, explicit gate as in the elements), `f1`/`f2` = `--gate-form 1|2` (p_h (phi^2 + f); the hyperbolic ANSYS cal5gp form). The seed pressure is ~1e-4 E_bio, so P >= 0.02 leaves g > 0.99; the law acts at P ~ 1e-4. `KEIO_PLAN.ja.md` §10 |
| files without `_v2` | development runs before 20 Oct evening: do not use |

Never put numbers of different rows of this table in one comparison without saying which is which.
