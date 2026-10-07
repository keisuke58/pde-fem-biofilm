# For the cloud session: review of the unattended ANSYS week (5 Oct 2026)

The user leaves for Paris on 6 Oct (about 09:30 CEST) for one week. IKMHIWI03
runs ANSYS on its own in that time; results arrive on this branch as JSON in
batches of ten. Your opinion is wanted **before 6 Oct 09:00**, while runs can
still be added or reordered on IKMHIWI03.

## What is queued (in this order)

1. tonight, already running: 16^3 case 3 (consumption 6), 24^3 one species
   beta = 0.02 (mesh series 8/16/24), then the rebuild with prop(35), a
   byte-identical regression check, and two prop(35) runs
   (`results/2026-10-05_ansys/`).
2. then 110 runs, `make_week_decks.py`, names in `results/2026-10-05_week/README.md`:
   - two species 8^3, cases 3 and 6, all with the Klempt 2024 stiffness
     (E = 10 Pa, nu = 0.49; the earlier two-species runs had the example 1000 MPa):
     consumption 0/2/4/6 x s 0.05/0.1/0.15/0.25/0.5/1.0, T* = 2 for
     s 0.05/0.15/0.5, phi_cap 0.85/0.95, start share 0.2/0.8, late start
     share (prop(35)) 0.2/0.5/0.8, c_ref 0.5/2, beta 0.01/0.05 (82 runs)
   - one species: beta 0.005/0.01/0.05/0.1 on 8^3 and 16^3, T* = 2, dt 0.0125
     on 16^3 (11 runs)
   - two species 16^3: consumption series for both cases, s 0.05/0.5, late
     0.2/0.8, beta 0.01/0.05, T* = 2 (15 runs)
   - one species 24^3: beta 0.01 and 0.05 (dt 0.0125) (2 runs)

Expected wall time about 1 to 1.5 days, so ANSYS is idle for the rest of the
week unless more is queued.

## Questions

1. Is anything missing that the thesis (ch. 4/5) or the 5 Oct open decisions
   (`COUPLING_STATUS.md`, "Open decisions") will need, and that cannot be run
   once the user is away? Anything in the list that is not worth the time?
2. What should fill the idle days? Candidates on IKMHIWI03's side:
   - two species on 24^3 (memory not yet checked; 32^3 one species did not fit
     in 31.5 GB),
   - a denser 16^3 grid (consumption x s) for the composition figures,
   - longer horizons (T* = 5) for one species, to see whether the seed stress
     keeps growing linearly.
3. For the s decision: which outputs would you use to choose s (T* = 2 runs,
   s x consumption grid)? If a specific quantity should be written per run
   (e.g. the seed-mean share over time from comp_trace), say so now: the JSON
   export currently keeps the last step of all_stress / nut_field and the
   element-220 stress history only; comp_trace stays on IKMHIWI03.
4. Is beta = 0.02 still the main value, given the beta series that will exist?

## Known limits of the queued set

- The first runs of the new decks start around 19:00 on 5 Oct; they are
  checked once (19:30) and again on 6 Oct morning before the user leaves.
- If the prop(35) rebuild fails, the old executable stays and the eight `late`
  runs do not use prop(35) (`_after18_1005.log` on IKMHIWI03 says which).
- A Windows Update reboot would stop the chain; none was pending on 5 Oct.
- Explicit phi update: dt beta / h^2 kept below about 0.1 (beta 0.1 on 16^3 and
  beta 0.05 on 24^3 use dt 0.0125); not checked beyond this estimate.

Answer by committing to this branch (a short file next to this one, or edits
to `make_week_decks.py` with the runs to add); IKMHIWI03 pulls it on 6 Oct
morning.
