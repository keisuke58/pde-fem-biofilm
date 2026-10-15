"""Step 1 of making the unified model the main tool: does growth as a nutrient reaction,
catalysed by the living cells, explain the measured volumes of Heine et al.?  (aggregates only,
data/heine_volume/; p3/unified_model.tex, Eq. growth)

In a packed biofilm the reaction rate is xidot = L0 * live(t) * (Y ln c - mu_b).  The nutrient
supply differs between the flow reactor (HOBIC, medium flows past) and static culture (medium
exchanged once a day), so the affinity A = Y ln c - mu_b is one number per cultivation, and only
the product k = L0 * A is identifiable from volumes.  live(t) is the measured community live
fraction (linear between days).  V = V0 (1 + xi)^3, fitted from Day 3 on (the drop from Day 1 to
Day 3 is medium exchange and detachment), weights 1/SD.

Tests:
 1. per cultivation, one k for the commensal and the dysbiotic community (V0 per condition);
    scheme (a) for comparison: the same with live = 1.
 2. prediction: k fitted on the commensal community only, then the dysbiotic volume predicted
    with that k (only V0 fitted, from its Day-3 value).
"""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data/heine_volume"


def load(fn):
    out = {}
    for r in csv.DictReader(open(D / fn)):
        out.setdefault((r["condition"], r["cultivation"]), {})[int(r["day"])] = r
    return out


VOL, LIVE = load("volume_summary.csv"), load("live_fraction_summary.csv")


def series(cond, cult, day0=3):
    v, l = VOL[(cond, cult)], LIVE[(cond, cult)]
    days = np.array(sorted(v), float)
    V = np.array([float(v[int(d)]["mean_biofilm"]) for d in days])
    S = np.array([float(v[int(d)]["sd_biofilm"]) for d in days])
    L = np.array([float(l[int(d)]["mean_biofilm"]) for d in days])
    m = days >= day0
    integ = np.array([np.trapezoid(L[(days >= day0) & (days <= d)], days[(days >= day0) & (days <= d)]) if d > day0 else 0.0
                      for d in days[m]])
    return days[m], V[m], np.maximum(S[m], 1e-3), integ, days[m] - day0


def fit(conds, cult, live=True, k_fixed=None):
    data = [series(c, cult) for c in conds]

    def model(k, v0, d):
        I = d[3] if live else d[4]
        return v0 * (1 + k * I) ** 3

    def res(p):
        k = p[0] if k_fixed is None else k_fixed
        v0s = p[1:] if k_fixed is None else p
        return np.concatenate([(model(k, v0, d) - d[1]) / d[2] for v0, d in zip(v0s, data)])
    p0 = ([0.01] if k_fixed is None else []) + [d[1][0] for d in data]
    r = least_squares(res, p0, bounds=(0, np.inf))
    return r.x, float(np.sum(r.fun ** 2)), sum(len(d[1]) for d in data), data


if __name__ == "__main__":
    for cult in ("HOBIC", "Static"):
        for live in (True, False):
            p, chi2, n, _ = fit(["Commensal", "Dysbiotic"], cult, live)
            print(f"{cult:6s} {'reaction, live' if live else 'scheme (a)    '}: k = {p[0]:.4f}/day, chi2 = {chi2:.2f} ({n} points, 3 parameters)")
        # prediction of the dysbiotic volume from the commensal fit
        for live in (True, False):
            p, chi2c, _, _ = fit(["Commensal"], cult, live)
            q, chi2d, nd, data = fit(["Dysbiotic"], cult, live, k_fixed=p[0])
            d = data[0]
            pred = q[0] * (1 + p[0] * (d[3] if live else d[4])) ** 3
            print(f"   predict dysbiotic with the commensal k ({'reaction' if live else '(a)     '}): chi2 = {chi2d:.2f} ({nd} points, V0 only); "
                  f"V(Day 21)/V(Day 3) measured {d[1][-1] / d[1][0]:.2f}, predicted {pred[-1] / pred[0]:.2f}")
