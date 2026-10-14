"""One-at-a-time sensitivity of the two-level model (unified_two_level_1d.py) to its new parameters.

Seed phi = 0.9 on x <= 0.25, outer liquid 1e-4, CH calibration, T* = 1, 60 steps.  Base:
lam = 1, k_alpha = 2, beta = beta_s = 0.02.  For each value: no antibiotic and antibiotic from
T* = 0.5 (alpha* = 0.65, b_i = 1, assumed).  Reported: alpha - 1 at the substrate, phi at x = 0,
0.3 and 1, living share at x = 0, and the antibiotic's effect on alpha - 1 and on phi(0)."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "JAXFEM", ROOT / "ansys_usermat/coupling"):
    sys.path.insert(0, str(p))
import unified_two_level_1d as t2  # noqa: E402
import unified_1d as u1  # noqa: E402

BASE = dict(lam=1.0, ka=2.0, beta=0.02, beta_s=0.02)
SWEEP = {"lam": [0.3, 1.0, 3.0], "ka": [0.5, 2.0, 8.0], "beta": [0.005, 0.02, 0.08], "beta_s": [0.005, 0.02, 0.08]}

if __name__ == "__main__":
    g0, th = u1.ch_case(); g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1
    x = np.linspace(0, 1, 41); phi_x = np.where(x <= 0.25, 0.9, 1e-4)
    print("param value | alpha-1(0) | phi(0) phi(0.3) phi(1) | living(0) | antibiotic: d(alpha-1) %, d phi(0)")
    for k, vals in SWEEP.items():
        for v in vals:
            if k != "lam" and v == BASE[k] and k != "ka":
                continue
            kw = dict(BASE, **{k: v}); res = {}
            for name, ast, on in (("none", 0.0, None), ("ab", 0.65, 0.5)):
                try:
                    outs, it, r, ns = t2.run(phi_x, g, th, ast=ast, ab_on=on, steps=60, **kw)
                    res[name] = outs[-1].reshape(41, t2.NV)
                except RuntimeError as e:
                    res[name] = None; print(k, v, name, "failed:", e)
            if res["none"] is None or res["ab"] is None:
                continue
            X, Y = res["none"], res["ab"]
            print(f"{k:6s} {v:6.3f} | {X[0, 13] - 1:.4f} | {X[0, 0]:.3f} {X[12, 0]:.3f} {X[-1, 0]:.1e} | "
                  f"{(X[0, 1:6] * X[0, 7:12]).sum():.3f} | {100 * ((Y[0, 13] - 1) / (X[0, 13] - 1) - 1):+.1f} %, {Y[0, 0] - X[0, 0]:+.3f}",
                  flush=True)
