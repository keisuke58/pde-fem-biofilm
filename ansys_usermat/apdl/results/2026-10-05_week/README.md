# Unattended ANSYS runs, from 5 Oct 2026 evening (one week, nobody at IKMHIWI03)

Decks: `ansys_usermat/apdl/make_week_decks.py` (run list and the reasons there).
Chain: `run_chain.ps1 -Name week_1005 ... -SkipDone -PushEvery 10`; the JSON files
arrive here in batches of ten, each batch one commit whose message carries the
chain log lines (rc per run). rc = 0 means ANSYS finished without errors; a run
with rc != 0 can still have a JSON (check its numbers before using them).

Read with `summarize_runs_json.py` (same format as `../2026-10-05_ansys/`).

Common to all runs: the partner's 2 mm cube, seed BIOFILM1, k_alpha = 1e-3 per T*
(field and Eq. 36), beta = 0.02 mm^2/T* unless the name says otherwise, dt = 0.025,
executable with the prop(35) fragment (5 Oct evening build).

## Names

Numbers drop the decimal point: `s015` = 0.15, `b005` = 0.05, `b0005` = 0.005,
`dt00125` = 0.0125, `T2` = T* 2.0, `s1` = 1.0.

Two species, `w<mesh>_<c3|c6>_g<consumption>_s<s>[...]`
- mesh 8 or 16 (8^3 / 16^3 elements), case 3 or case 6 of Klempt et al. 2026
- `g` = CONSUMPTION11 (0 = no consumption: c = 1 everywhere). `g1` is Klempt
  2024 Table 2 scaled to the 2 mm cube like beta (g/d with d = MY_DIFF1 = 1;
  see make_week_decks.py); 2/4/6 are example inputs
- `s` = clock of the point model (prop(31)), an assumption
- optional: `_cap` phi_cap (prop(32), default 0.9), `_chi` start share (prop(30),
  default 0.5), `_late` start share of points reached later (prop(35), default
  off = prop(30)), `_cref` c_ref (prop(33), default 1), `_b` beta, `_T` end time
  (default T* = 1)
- stiffness: Klempt 2024 (E = 10 Pa, nu = 0.49), unlike the two-species runs in
  `../2026-10-05_ansys/` (example input 1000 MPa); composition does not depend on it

One species, `w<mesh>_1sp_b<beta>[_dt<dt>][_T<T>]`
- E = 10 Pa, nu = 0.49, end time T* = 1.1 (unless `_T`)
- beta = 0.1 on 16^3 and 0.05 on 24^3 use dt 0.0125 (explicit phi update, keeps
  dt beta / h^2 below about 0.1)

One species variants for ch4 (cloud review, `../../CLOUD_REVIEW_1005_ANSWERS.md`):
`w<mesh>_1sppartner_b002` the partner's own growth variable (prop(28) = 3),
`w<mesh>_1spnu0499_b002` nu = 0.499, `w<mesh>_1spEex_b002` the example stiffness
(1000 MPa, nu = 0.3, linear blend).

`share_history` in each two-species JSON: seed and layer 1/2 mean, min, max of
phi_1/(phi_1+phi_2), mean amount, number of points, per substep (from
comp_trace). The executable of 5 Oct evening traces every point with
phi >= 0.01; `points` in the JSON says which rule a run had. Nut1 over time is
not traced (only the last step, in nut_field).

## Order (137 runs)

0. one species, ch4 variants, 8^3 and 16^3 (6 runs)
1. two species 8^3 (90 runs, about 3-15 min each)
2. one species 8^3 / 16^3 (11 runs)
3. two species 16^3 (17 runs, about 30-60 min each)
4. one species 24^3, beta 0.01 and 0.05 (2 runs, 2-4 h each)
5. idle days: two species 24^3 case 6, consumption 1, 6 and 0 (memory unchecked;
   a FATAL stop is not retried), one species T* = 5 on 8^3 / 16^3, six more
   16^3 two-species runs
