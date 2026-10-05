# Answers to CLOUD_REVIEW_1005_WEEK.md (cloud session, 5 Oct evening)

Based on the JSON in `results/2026-10-05_ansys/` and the map of the ch4 values
to replace (`BETA002_THESIS_MAP.md`). Run names below are suggestions; adapt
them to `make_week_decks.py`.

## 1. Missing / not worth the time

Add (all short, 8^3 unless stated; put them early, they feed ch4 directly):

| # | run | why |
|---|---|---|
| A | one species, **the partner's own growth variable**, beta = 0.02 (same deck as the 2 Oct `ds_kl_partner`-type run, only `MY_BETA1=0.02`) | ch4 figure of the seeded element over time and "The whole model" (ratio 0.455) compare Eq. 36 with the partner's variable; today only beta = 1e-4 exists for the partner's variable |
| B | one species, beta = 0.02, **nu = 0.499** (`POISSON_BIO=0.499`) | the incompressibility statement in ch4 ("mean stress changes by < 10 % between 0.49 and 0.499") needs the beta = 0.02 value; nu = 0.45 exists |
| C | one species, beta = 0.02, **framework example stiffness** (E = 1000 MPa, nu = 0.3, linear blend, as the 2 Oct comparison) | optional: the sentence comparing the published and the example stiffness (factor 6e-9, cosine 0.94) is at beta = 1e-4 only |

Not worth the time (drop if the queue gets long):
- `cref 0.5 / 2`: only c / c_ref enters, so c_ref = 2 repeats the consumption
  series (half the nutrient) and c_ref = 0.5 only clips at 1. Four runs.
- `s = 1.0` with consumption 0 and 2 on 8^3: at s = 0.5 the seed share is
  already 0.07 (case 6), i.e. the outcome is reached; s = 1.0 adds nothing
  except at consumption 4/6 (keep those two).

No re-run needed: the exact-solution check (sinh/cosh) is a check of the time
stepping in a gradient-free region and stays at beta = 1e-4, stated as such.

## 2. Idle days, in this order

1. **Two species on 24^3, case 6, consumption 6, beta = 0.02** (and its
   no-nutrient control). The composition range in the seed widens from 8^3 to
   16^3 (0.117-0.459 to 0.070-0.466) because the nutrient face is resolved;
   a third mesh is what the thesis needs to call it converged. Check memory
   first with a short run (T* = 0.1): one species 24^3 needed about 7 GB; the
   point model adds the Python server, not PARDISO memory.
2. **T* = 5, one species, beta = 0.02, 8^3 and 16^3**: whether the seed stress
   keeps growing linearly and how far phi spreads (the "Magnitude" limit in
   ch4). Check dt beta / h^2 as for the other runs.
3. Only then the denser 16^3 consumption x s grid.

## 3. Outputs for the s decision

The open decision note says the share at T* = 1 follows s T*. ANSYS output
alone cannot fix s: s links the point model's clock to the field's T*, and
choosing it needs a time scale from data (the TMCMC calibration window of the
point model is 0.25 model time = 21 days of the Heine experiment; the field's
T* has no stated physical time in Klempt 2024). What the thesis needs from the
runs is the shape of the dependence, so please export per run:

- **the seed-mean share phi_1/(phi_1+phi_2) over time** (and seed min/max, and
  the means of layers 1 and 2), one row per coupling step, from comp_trace or
  nut_field per substep. A few kB per run, enough to plot share against s T*
  and to test that s and T* enter only as s T* (the T* = 2 runs with s = 0.05,
  0.15, 0.5 against T* = 1 with s = 0.1, 0.25 (nearest), 1.0);
- the seed-mean amount phi_1 + phi_2 and Nut1 over time (same rows).

With these I can write s as an assumption with one figure: share against s T*
for consumption 0 and 6, both cases.

## 4. beta = 0.02 as the main value

Yes. It is the paper's value converted to this geometry (stated as an
assumption), and it is the only value tested so far for which the seed mean
stress (0.6 % between 8^3 and 16^3) and the first-layer tension (0.7 %) are
mesh converged. The queued beta series (0.005 to 0.1) then goes into ch4 as
one sensitivity line or small table: seed mean stress and first-layer stress
against beta, 8^3 and 16^3. If the seed stress turns out to scale strongly
with beta, the thesis says so; the choice of 0.02 does not depend on it.

Note for the thesis (already in BETA002_THESIS_MAP.md): with beta = 0.02 the
composition outside the seed depends on the start share of newly reached
points (case 6 has two stable outcomes, split between 0.5 and 0.7), so the
prop(35) runs matter; keep all three (0.2/0.5/0.8).
