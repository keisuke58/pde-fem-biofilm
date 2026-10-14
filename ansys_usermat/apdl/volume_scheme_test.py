"""Do the measured volumes of Heine et al. follow growth from phi (scheme a) or from the
living biomass (scheme c)?  Aggregates only (data/heine_volume/).

Packed biofilm (phi ~ 1): alpha(t) = 1 + k * int_{d0}^{t} g dt,  V = V0 * alpha^3,
  (a) g = 1,   (c) g = L(t), the measured community live fraction (linear in time between days).
Fitted on the growing part only (from the day of the smallest volume on), weights 1/SD.
Per condition: V0 and k free (2 parameters in both schemes).  Shared: one k for all
conditions of a cultivation, V0 per condition.  The current growth law cannot shrink,
so the falling parts are left out; that choice is stated with the result."""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data/heine_volume"


def load(fn, key):
    out = {}
    for r in csv.DictReader(open(D / fn)):
        out.setdefault((r["condition"], r["cultivation"]), {})[int(r["day"])] = r
    return out


VOL = load("volume_summary.csv", None)
LIVE = load("live_fraction_summary.csv", None)


def series(cond, cult):
    v, l = VOL[(cond, cult)], LIVE[(cond, cult)]
    days = sorted(v)
    V = np.array([float(v[d]["mean_biofilm"]) for d in days])
    S = np.array([float(v[d]["sd_biofilm"]) for d in days])
    L = np.array([float(l[d]["mean_biofilm"]) for d in days])
    return np.array(days, float), V, S, L


def growth_integral(days, L, d0, scheme):
    """int_{d0}^{d} g dt at each day >= d0, g linear between measured days."""
    g = np.ones_like(L) if scheme == "a" else L
    out = []
    for d in days:
        if d < d0:
            out.append(np.nan); continue
        m = (days >= d0) & (days <= d)
        x, y = days[m], g[m]
        out.append(np.trapezoid(y, x) if len(x) > 1 else 0.0)
    return np.array(out)


def segment(cond, cult):
    days, V, S, L = series(cond, cult)
    i0 = int(np.argmin(V[:-1])) if np.argmin(V) == len(V) - 1 else int(np.argmin(V))
    m = days >= days[i0]
    return days, V, S, L, days[i0], m


def fit_one(cond, cult, scheme):
    days, V, S, L, d0, m = segment(cond, cult)
    I = growth_integral(days, L, d0, scheme)[m]
    Vm, Sm = V[m], np.maximum(S[m], 1e-3)
    res = lambda p: (p[0] * (1 + p[1] * I) ** 3 - Vm) / Sm
    r = least_squares(res, [Vm[0], 0.01], bounds=([1e-6, 0], [10, 10]))
    return r.x, float(np.sum(r.fun ** 2)), int(m.sum())


def fit_shared(cult, scheme, conds):
    segs = [segment(c, cult) for c in conds]
    Is = [growth_integral(s[0], s[3], s[4], scheme)[s[5]] for s in segs]

    def res(p):
        k = p[0]
        return np.concatenate([(p[1 + j] * (1 + k * Is[j]) ** 3 - s[1][s[5]]) / np.maximum(s[2][s[5]], 1e-3)
                               for j, s in enumerate(segs)])
    p0 = [0.01] + [s[1][s[5]][0] for s in segs]
    r = least_squares(res, p0, bounds=([0] + [1e-6] * len(segs), [10] + [10] * len(segs)))
    return r.x, float(np.sum(r.fun ** 2)), sum(int(s[5].sum()) for s in segs)


if __name__ == "__main__":
    conds = ["Commensal", "Dysbiotic"]
    for cult in ["HOBIC", "Static"]:
        for cond in conds:
            days, V, S, L, d0, m = segment(cond, cult)
            line = f"{cond[:3]}-{cult[:3]} from day {int(d0)} ({int(m.sum())} points): "
            for sc in "ac":
                p, chi2, n = fit_one(cond, cult, sc)
                line += f" ({sc}) k={p[1]:.4f}/day chi2={chi2:.2f}"
            print(line)
        for sc in "ac":
            p, chi2, n = fit_shared(cult, sc, conds)
            print(f"  shared k, {cult}, ({sc}): k={p[0]:.4f}/day chi2={chi2:.2f} ({n} points, {1+len(conds)} parameters)")
