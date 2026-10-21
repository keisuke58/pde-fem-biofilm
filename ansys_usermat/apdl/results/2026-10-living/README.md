# ANSYS runs of the coupling agreed with Oliver (20 Oct 2026): which results to use

Generations, oldest first. **Use only the current generation in the thesis, the slides and mails.** Never put
numbers of different generations in one table without saying which is which.

| prefix | what | status |
|---|---|---|
| `lv8_`, `trace_*` (A-19) | first runs | **do not use**: L frozen at its first value (A-22) |
| `v2_`, `v2s_`, `sens_`, `i8*` | v2 set-up | **do not use**: L frozen (A-22) |
| `v2f_`, `fsens_`, `fthr_`, `sym_*`, `dbg_*` | L fixed | **do not use for results**: phi advanced 10x too fast (inner USSFin step 0.0025, A-25). `sym_*` stay valid only as the neighbour/symmetry test (A-24) |
| `u1_` (A-27..A-31) | one inner step (`--dt-ussfin 0.1`), quasi-static c and a | **main ANSYS results** (21 Oct, later: with d of Klempt 2024 Table 2, 1e8 mm^2/T* on the 2 mm cube, the time derivative changes the results by < 0.005, Python) |
| `u2_` (addendum 23, optional) | `NUT_TRANSIENT = 1` | only as a check with D = 1e8 (should equal u1_), or for a slow-diffusion sensitivity |

Figures made from these runs: `assets/anim_u1_16_*.gif`, `assets/frames_u1_16_*.pdf` (quasi-static u1_, comparison only;
remake from u2_), `assets/fig_coupled_examples.pdf` (Python, quasi-static, before the start rule: do not use).
`assets/fig_contact_zone.pdf` is from the time-derivative Python runs (`_transc1`): current.
