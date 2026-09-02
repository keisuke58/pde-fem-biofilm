# Finding — the φ→E bridge saturates over much of the healthy composition space

Found 2026-09-02, while looking for failures that are silent. **Not changed**:
the clamp is deliberate and may be the right physical bound. What is worth
knowing is how often it binds, because when it does, two different conditions
report the same stiffness and nothing says so.

## The mechanism

`material_models.compute_E_composite` ends with

```python
    # Clamp to physical range
    return np.clip(E_eff, E_MIN_PA, E_MAX_PA)      # 10 Pa, 1000 Pa
```

and it is the production bridge — `ansys_usermat/coupling/composition_to_material.py`
line 110 calls it to turn a composition into the `(C10, C01, D1)` a deck
carries. (Lineage 2's `compute_E_di` is a different bridge; see
`RESEARCH_MODEL.md` §6, and do not mix them in the write-up.)

At the current calibration (`MECH_E0 = 221340 Pa`, `α = 2.07`) the mechanistic
model produces values well above the cap for EPS-rich compositions:

| composition | model output | reported | over cap |
|---|---|---|---|
| all commensal (*S. oralis*) | 3127 Pa | **1000 Pa** | 3.1× |
| even five-species mix | 2991 Pa | **1000 Pa** | 3.0× |
| *S. oralis* + *A. naeslundii* | 7615 Pa | **1000 Pa** | 7.6× |
| mostly *P. gingivalis* | 256 Pa | 256 Pa | — |
| all *P. gingivalis* | 0 Pa | 10 Pa | — |

## Why it matters, and what it is not

The headline of this work is a **comparison between conditions**. Two
compositions whose modelled stiffness differs by 2.4× — the even mix against
*S. oralis* + *A. naeslundii* — both come out at exactly 1000 Pa, so the stress
fields they produce are identical and the comparison shows no difference where
the model says there is a large one. No warning is emitted; a clip is silent
by construction.

This is **not** an argument that the clamp is wrong. `E_MAX_PA` is commented as
the commensal upper limit, i.e. an intended physical bound, and a mechanistic
power law with an exponent above 2 will run away without one. The finding is
narrower and entirely about disclosure: at this calibration the bound is
*active over much of the healthy half of composition space*, not a rare guard
against outliers.

## What is not settled

**Do the actual CLSM conditions saturate?** The compositions above are
synthetic corners. Whether the conditions this thesis reports sit above the cap
cannot be determined here — the extract CSVs are not in this checkout (the same
files `JAXFEM/audit_all.py` skips on). This is the one question that decides
whether the finding affects the reported numbers at all, and it is a one-command
check once the data is at hand:

```
python -c "
import numpy as np, material_models as mm
E = mm.compute_E_composite(phi)          # phi = your conditions, (n, 5)
print((np.isclose(E, mm.E_MAX_PA) | np.isclose(E, mm.E_MIN_PA)).sum(),
      'of', len(E), 'conditions are on a bound')
"
```

`tests/test_e_saturation.py` pins the behaviour so it cannot change unnoticed.

## If it turns out to bind for the reported conditions

Three options, in increasing order of work — a decision, not a
recommendation:

1. **Report it as a limitation.** Cheapest, and honest: state that the bridge
   saturates and that differences among healthy-side conditions are therefore
   lower bounds.
2. **Raise `E_MAX_PA`** to whatever the literature supports for commensal
   biofilm. Changes every reported number; needs the citation first.
3. **Recalibrate `MECH_E0` / `α`** so the model lands inside the bound over the
   observed range. Most defensible, most work, and it re-opens the Pattem 2018
   calibration.
