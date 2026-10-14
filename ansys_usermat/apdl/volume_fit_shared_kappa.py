"""Scheme (c) with one kappa for all conditions, and with decay of the dead biomass,
fitted to the measured volumes of Heine et al. (aggregates, data/heine_volume/).

Growth at a seed Gauss point of the element runs (point model as in scheme b, its own
phi_i, psi_i; results/2026-10-cal5/w8_<C>_cal_k0882_pm.json):
    d ln(alpha)/dT* = kappa * chi * sum_i phi_i psi_i  -  k_d * chi * sum_i phi_i (1 - psi_i)
    V = V0 * alpha^3,   Day = 1 + 20 T*  (T* = 1 is Day 21, the element's time link)
chi = phi_3D (phi_used).  The point model does not depend on alpha, so the traces hold for
any kappa and k_d.  kappa and k_d shared by the conditions, V0 per condition; the decay
term is an assumption (dead cells broken down at rate k_d), not from a paper.
Only the conditions with a MAP in the element (CH, CS; flow reactor) are used."""
import csv, json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "ansys_usermat/apdl/results/2026-10-cal5"
COND = {"CH": ("Commensal", "HOBIC"), "CS": ("Commensal", "Static")}


def trace(c):
    pm = json.load(open(R / f"w8_{c}_cal_k0882_pm.json")); C = {k: i for i, k in enumerate(pm["columns"])}
    rows = pm["seed_trace"]; rows = [r for r in rows if r[C["ip"]] == rows[0][C["ip"]] and r[C["elem"]] == rows[0][C["elem"]]]
    d = {int(r[C["isubst"]]): r for r in rows}; k = sorted(d); rows = [d[i] for i in k]
    t = np.array([0.0] + list(np.array(k) / max(k)))
    phi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(1, 6)] for r in rows])
    psi = np.array([[float(r[C[f"gnew_{i}"]]) for i in range(7, 12)] for r in rows])
    chi = np.array([float(r[C["phi_used"]]) for r in rows])
    live = chi * (phi * psi).sum(1); dead = chi * (phi * (1 - psi)).sum(1)
    pad = lambda x: np.concatenate([[x[0]], x])
    return t, pad(live), pad(dead)


def cumint(t, g):
    return np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])


def volumes(cond, cult):
    rows = [r for r in csv.DictReader(open(ROOT / "data/heine_volume/volume_summary.csv"))
            if r["condition"] == cond and r["cultivation"] == cult]
    d = np.array([float(r["day"]) for r in rows]); V = np.array([float(r["mean_biofilm"]) for r in rows])
    S = np.array([float(r["sd_biofilm"]) for r in rows])
    return d, V, S


import sys
DAY0 = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0   # Day 1 -> 3 is a drop in all conditions (medium exchange, detachment)
DATA = {}
for c, (cond, cult) in COND.items():
    t, live, dead = trace(c)
    d, V, S = volumes(cond, cult)
    m = d >= DAY0; d, V, S = d[m], V[m], S[m]
    Td = (d - 1) / 20
    IL = np.interp(Td, t, cumint(t, live)); ID = np.interp(Td, t, cumint(t, dead))
    DATA[c] = (IL - IL[0], ID - ID[0], V, S)


def model(p, c, decay):
    kappa, kd = p[0], (p[1] if decay else 0.0)
    v0 = p[2 + list(COND).index(c)] if decay else p[1 + list(COND).index(c)]
    IL, ID, _, _ = DATA[c]
    return v0 * np.exp(3 * (kappa * IL - kd * ID))


def fit(decay):
    def res(p):
        return np.concatenate([(model(p, c, decay) - DATA[c][2]) / DATA[c][3] for c in COND])
    p0 = [0.5, 0.5] + [DATA[c][2][0] for c in COND] if decay else [0.5] + [DATA[c][2][0] for c in COND]
    lo = [0, 0] + [1e-3] * len(COND) if decay else [0] + [1e-3] * len(COND)
    r = least_squares(res, p0, bounds=(lo, [50] * len(p0)))
    return r.x, float(np.sum(r.fun ** 2))


if __name__ == "__main__":
    for decay in (False, True):
        p, chi2 = fit(decay)
        n = sum(len(DATA[c][2]) for c in COND)
        print(f"decay={decay}: params {np.round(p, 4)}  chi2 {chi2:.2f}  ({n} points, {len(p)} parameters)")
        for c in COND:
            print("   ", c, "measured", np.round(DATA[c][2], 3), "model", np.round(model(p, c, decay), 3))
