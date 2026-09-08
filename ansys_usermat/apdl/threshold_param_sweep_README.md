# Critical-coupling threshold: dependence on dt / Penalty1 / MaxGrowth

2026-09-08 (last real-ANSYS day). Follow-up to `N3_GROWTH_TRIAL.md`'s
bisected critical coupling threshold for the species-3↔4 interaction
(baseline `MaxGrowth=100/dt=0.1/Penalty1=5/HalfVelo=0.1`: stable at
`0.0005`, diverges at `0.0007` — narrowed here to diverges at `0.0006`
too, so the true baseline threshold sits in `(0.0005, 0.0006)`). That
file's own note flagged this as the natural next step, not attempted that
session: "test coupling magnitudes below each new threshold, not just
retry the same 0.05 value" for a 10x change to each of `dt`, `Penalty1`,
`MaxGrowth`.

Reused the three existing parameter-variant decks from that session
(`ds_oliver_wired_n5trial_inter34only_{dt01,pen50,lowgrowth}.dat`, each a
10x change to one constant, previously only tested at the
far-above-threshold coupling `0.05`) and re-targeted their
`INTERACTION34`/`INTERACTION43` at values near the baseline threshold
instead, to find each variant's own threshold. All runs: real ANSYS,
`F:\biofilm_upf_wired`, `-smp -np 1`, unmodified `Ussfin`/`Usermat` build
from 2026-09-07 (no recompile needed). Full PASS/FAIL lines: `threshold_param_sweep_SUMMARY.txt` in this same
directory; the `ds_oliver_wired_n5trial_inter34only_*_at0p*.dat` deck
files sit alongside it (raw ANSYS output logs were NOT kept — 20-27 MB
each per run, impractical to commit — the summary file's grepped lines
are the record).

## Results

| Variant | Change from baseline | Coupling tested | Result |
|---|---|---|---|
| baseline | — | 0.0006 | **diverges** (new data point, not tested in the original bisection) |
| dt01 | `dt`: 0.1 → 0.01 (10×smaller) | 0.0006 | diverges |
| dt01 | same | 0.0005 (baseline's known-*pass* point) | **diverges** |
| pen50 | `Penalty1`: 5 → 50 (10×larger) | 0.0006 | diverges |
| pen50 | same | 0.0005 | passes (0 errors) |
| lowgrowth | `MaxGrowth13/14`: 100 → 10 (10×smaller, species 3/4 only) | 0.0006 | passes |
| lowgrowth | same | 0.005 (10× the baseline threshold) | **passes** |

## Reading it

- **`MaxGrowth` is the dominant lever, and it's a large one.** Cutting the
  coupled species' own growth rate 10× raises the stable coupling ceiling
  by *at least* ~10× (0.0006 → ≥0.005, not yet bisected further — only a
  lower bound on how much higher it goes). Physically sensible: a species
  growing less aggressively on its own has more headroom before an added
  interaction term pushes it over the edge — consistent with
  `N3_GROWTH_TRIAL.md`'s root-cause finding that species 3's own
  *uncoupled* trajectory is already close to unstable under the baseline
  constants.
- **`Penalty1` (10× larger) does not move the threshold much**, if at
  all — still fails at 0.0006 and still passes at 0.0005, the same
  bracket as the unmodified baseline. Not the lever to reach for.
- **`dt` (10× smaller) makes it *worse*, not better** — diverges even at
  0.0005, where baseline (and pen50) are stable. This is *not* read as
  "smaller timesteps destabilize this scheme" in a generic numerical
  sense: `TIME` is held fixed at `1.1` in both baseline and `dt01`
  (`deltim` 0.1 → 11 substeps; 0.01 → 110 substeps), so `dt01` is not the
  same physical duration taken in finer numerical steps — it is **the
  same coupling strength applied for 10× as many explicit growth-update
  increments**, which more directly explains a lower apparent threshold
  than any claim about numerical stiffness. Worth remembering when this
  is picked back up: to isolate the numerical-dt-vs-total-duration
  question, `TIME` would need to shrink along with `dt` (i.e. compare
  equal numbers of substeps, not equal elapsed time), which was not done
  here.

## Not done

- No further bisection to pin `lowgrowth`'s actual new threshold (only
  bracketed as `[0.0006, 0.005]` from below/above at the two points
  tested).
- No combined changes (e.g. lower `MaxGrowth` *and* higher `Penalty1`
  together) — `N3_GROWTH_TRIAL.md` already noted this as the natural
  extension once single-parameter effects are known, still true here.
- The `dt`-vs-`TIME` confound above is flagged, not resolved.
