# results/2026-10-living/archive (22 Oct 2026)

ANSYS runs of the living-share coupling that are **not to be used** for results, kept for the record. Why each
generation is superseded: `../README.md` (generations table) and `../../../CONFIRMED_2026-10-22.ja.md` §3.

| prefix | why |
|---|---|
| `lv8_`, `tsens_`, `v2_`, `v2s_`, `sens_`, `i8*` | the living share L was read only once (frozen), A-22 |
| `v2f_`, `fsens_`, `fthr_`, `v2t_`, `dbg_`, `ussfin_L_trace_dbg_*` | phi advanced ten times too fast (inner USSFin step), A-25 |

The diagnostic scripts that read them (`check_living_fields.py`, `summarize_living.py`, `fig_living_trace_check.py`,
`nem_neighbour_symmetry.py`) look here. Kept in `..`: the test runs `sym_*` (valid as the neighbour/symmetry test) and every
run of the generations listed there as checks or main results.
