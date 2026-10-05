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
- `g` = CONSUMPTION11 (0 = no consumption: c = 1 everywhere), example input
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

## Order

1. two species 8^3 (82 runs, about 3-15 min each)
2. one species 8^3 / 16^3 (11 runs)
3. two species 16^3 (15 runs, about 30-60 min each)
4. one species 24^3, beta 0.01 and 0.05 (2 runs, 2-4 h each)
