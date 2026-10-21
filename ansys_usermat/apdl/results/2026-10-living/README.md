# ANSYS runs of the coupling agreed with Oliver (20 Oct 2026): which results to use

Generations, oldest first. **Use only the current generation in the thesis, the slides and mails.** Never put
numbers of different generations in one table without saying which is which.

| prefix | what | status |
|---|---|---|
| `lv8_`, `trace_*` (A-19) | first runs | **do not use**: L frozen at its first value (A-22) |
| `v2_`, `v2s_`, `sens_`, `i8*` | v2 set-up | **do not use**: L frozen (A-22) |
| `v2f_`, `fsens_`, `fthr_`, `sym_*`, `dbg_*` | L fixed | **do not use for results**: phi advanced 10x too fast (inner USSFin step 0.0025, A-25). `sym_*` stay valid only as the neighbour/symmetry test (A-24) |
| `u1_` (A-27..A-31) | one inner step (`--dt-ussfin 0.1`), **quasi-static** c and a (0 = d lap c - ...) | valid, but only as the **quasi-static comparison** |
| `u2_` (addendum 23, after the trip) | `NUT_TRANSIENT = 1`: dc/dt and da/dt as in Eq. 35 and Oliver's sketch | **the main ANSYS results** once they exist |
