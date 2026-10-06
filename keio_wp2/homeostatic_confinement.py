#!/usr/bin/env python3
"""Keio WP2, step 3: the homeostatic-pressure growth law needs confinement.

In homeostatic.py the seed sits in a void 1000 times softer than the biofilm
(f = 1e-3, the floor of the runs), expands almost freely and builds only
about 5e-4 Pa. Measured colonies stop growing at pressures of about 10 kPa
when they are confined (Chu et al. 2018, CHU2018_NOTES.md), e.g. in a gel
(Vibrio cholerae in agarose) or against a substrate. Here the surrounding
stiffness is raised, E_out = f E, f = 1e-3 ... 1, and the homeostatic pressure
is fixed relative to the fully confined pressure of the biofilm,

    p_h = P_h * K (k_alpha T*),     K = E / (3 (1 - 2 nu)),

so that the result depends only on f and the dimensionless P_h. Growth law

    alpha_dot = k_alpha phi max(0, 1 - p / p_h),

phi as in Klempt 2024 Eq. 34 without front term and nutrient. 3D octant,
spherical seed (stress_to_growth_3d.Mesh), 16^3, dt = 0.02, T* = 1.

    python keio_wp2/homeostatic_confinement.py -> keio_wp2/results_confinement.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import homeostatic as H  # noqa: E402
import stress_to_growth_3d as S  # noqa: E402

K_BIO = S.E0 / (3 * (1 - 2 * S.NU))
P_H = 0.25                                       # p_h / (K k_alpha T*)


def run(f, p_h, n=16, dt=0.02, t_end=1.0):
    old = S.FLOOR
    S.FLOOR = f
    try:
        M = S.Mesh(n, "sphere")
        phi = M.inside.astype(float)
        alpha = np.ones_like(phi)
        hist = []
        for k in range(int(round(t_end / dt))):
            p = H.mechanics_p(M, phi, alpha - 1)
            g = np.ones_like(p) if p_h is None else np.clip(1 - p / p_h, 0.0, 1.0)
            phi = np.clip(phi + dt * (S.BETA * M.lap(phi) + S.K_ALPHA * alpha), 0, 1)
            alpha = alpha + dt * S.K_ALPHA * phi * g
            ins = M.inside
            hist.append(((k + 1) * dt, float((alpha - 1)[ins].mean()), float(p[ins].mean()),
                         float(p[ins].max()), float(g[ins].mean())))
    finally:
        S.FLOOR = old
    return dict(f=f, p_h=p_h, hist=hist)


def main():
    t0 = time.time()
    p_h = P_H * K_BIO * S.K_ALPHA
    out = {"P_h": P_H, "p_h": p_h, "K_bio": K_BIO, "runs": []}
    for f in (1e-3, 1e-2, 1e-1, 1.0):
        for ph in (None, p_h):
            r = run(f, ph)
            out["runs"].append(r)
            t, d, pm, px, g = r["hist"][-1]
            print(f"f={f:7.0e} p_h={'none' if ph is None else f'{ph:.3e}'}: alpha-1 {d:.3e}, p mean {pm:.3e} "
                  f"max {px:.3e} Pa (p_h {p_h:.3e}), g {g:.2f} [{time.time() - t0:.0f} s]", flush=True)
            (HERE / "results_confinement.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
