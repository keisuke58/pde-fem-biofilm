# ANSYS runs of the coupling agreed with Oliver (20 Oct 2026): which results to use

Generations, oldest first. **Use only the current generation (`u8_`) in the thesis, the slides and mails.** Never put
numbers of different generations in one table without saying which is which. The confirmed facts behind this table:
`../../CONFIRMED_2026-10-22.ja.md` (updated 22 Oct 2026).

| prefix | build (`F:\`) | what | status |
|---|---|---|---|
| `lv8_`, `trace_*`, `tsens_` (A-19) | living | first runs | **do not use**: L frozen at its first value (A-22); moved to `archive/` |
| `v2_`, `v2s_`, `sens_`, `i8*` | living | v2 set-up | **do not use**: L frozen (A-22); in `archive/` |
| `v2f_`, `fsens_`, `fthr_`, `v2t_`, `dbg_*` | living4 | L fixed | **do not use**: phi advanced 10x too fast (inner USSFin step 0.0025, A-25); in `archive/` |
| `sym_*`, `symt_*` | living4 | neighbour / symmetry test | test only (A-24), not results |
| `u1_` (A-27..A-31) | living4 | one inner step, quasi-static c, values before the paper-value decision | **superseded** |
| `u3_` (A-32) | living4 | paper front r, no stabilisation | **do not use**: elements distort at 8^3 and 16^3 |
| `u3a_`, `u3a5_`, `u3t_` (A-33..A-39) | living4/5 | artificial diffusion r h/4 at 8^3 (quasi-static / dc/dt), incl. the time-step series | 8^3 checks only; phi leaks through the walls |
| `u4_` (A-35 16^3, A-36 24^3) | living5 | production set-up of 21 Oct (dc/dt, d = 1e8, --artdiff, increment 0.00125) | **superseded**: phi leaks through the walls (surface points, A-39). `u4_16_*` stay only as the surface-point sensitivity; `u4_24_*` stopped |
| `u6_`, `u6t_`, `u6d_` (A-38, A-39) | living6 | NUT1_INIT_ALL (c(0) = 1 everywhere), front dump | 8^3 checks; same as `u3a_` (leak) |
| `u6s_`, `u5_` (A-39, A-40) | living6 | + SURF_NEAREST (leak / 4) | intermediate; `u5_16_sep_ab` is addendum 30's u4s_16_sep_ab |
| `u7c_`, `u7n_`, `u7p*` (A-42, A-44) | living7 | + LAP_FILE (conservative Laplacian) at 8^3; `u7n_` = short NEM-dump runs, `u7p*` = timing tests (no results) | superseded: neighbours by element order (asymmetric inside at 16^3 and 24^3, A-44); `u7_` was never run |
| **`u8_`** (A-45) | **living8** | living7 + NEM neighbours by distance; SURF_NEAREST + LAP_FILE, dc/dt d = 1e8, --artdiff N, increment 0.00125 (`u8c_` = 8^3 check, `u8n_` = dumps) | **current main ANSYS results** (front term stabilised by r h/4: stresses low, A-57) |
| `lw16_`, `lw24_` (A-59, A-62) | living10 | u8 + Lax-Wendroff front + sensor-switched added diffusion (r h/8) | 16^3 = second ANSYS result for the stress; 24^3 unstable (corners) |
| `up8_`, `up16_` (A-64) | living11 | u8 + upwind MLS gradients in the front term | stable on 16^3; 24^3 (`up24_`) running |
| `rt*`, `mr*`, `mt*`, `mb16`, `mh16`, `mms*b` (A-65, A-70) | living12-14 | manufactured solution (MMS_SRC): raw element, mirrored neighbours (`mt`/`mb`/`mh`: beta* 0.2 / 0.1 / half time step) | verification only |
| `mu8_` (A-72) | living14 | mirrored neighbours + upwind in the real example | 8^3 only; 16^3 unstable at a corner (not exported) |
| `gci_a52.*` | | grid convergence index of all codes | current |

Moved to `archive/` (24 Oct 2026): every generation before `u8_` (`u1_` to `u7*`, `sym*`), the Klempt 2024 Sect. 4.2
attempts (`k24*`, not reproduced), the face-constraint tests of A-69 (`f*m*`), the wall-sensor tests (`lw8p0`,
`lw8w*`, A-63) and `x10a`. `checks_a38_a39/summary_a38.py` and the default runs of `anim_living_fields.py` read
some of them: pass the `archive/` path.

Deck options of `u8_`: `--front 1 --v2 --start-step 0.0025 --dt-ussfin 0.1 --cons 1 --cons-ab 1 --phi-min 0.001
--artdiff N --deltim 0.00125 --transient-c 1e8 --nut1-init-all --surf-nearest --lap-file`, with `lap_cons.txt` of the
same mesh in the work directory (`../../checks_a38_a39/build_lap_cons_a42.py` from the `nem_dump.txt` of a `--nem-dump` run).

Figures made from these runs: `assets/anim_u1_16_*.gif`, `assets/frames_u1_16_*.pdf` (u1_, superseded: remake from u7_),
`assets/fig_coupled_examples.pdf` (Python, quasi-static, before the start rule: do not use).
`assets/fig_contact_zone.pdf` is from the time-derivative Python runs (`_transc1`).
