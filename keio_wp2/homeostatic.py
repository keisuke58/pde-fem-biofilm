#!/usr/bin/env python3
"""Keio WP2, step 2: the two ways forward after the first prototype.

(A) The Eq. 30 term at large Pi, as a point model (0D) for the seed surface:

        phi_dot = k_alpha alpha - mu_star c phi delta^2,   delta_dot = k_alpha phi,   delta = alpha - 1.

    For Pi >> 1 phi relaxes to phi_eq = k_alpha alpha / (mu_star c delta^2), and
    delta_dot = k_alpha phi_eq gives

        delta(t) ~ (3 k_alpha^2 t / (mu_star c))^(1/3)      (growth slows as t^(1/3)),

    but only for t >> 1/k_alpha. On T* = 1 the first transient sets the result:
    phi is removed once mu_star c delta^2 t ~ 1 with delta ~ k_alpha t, which gives
    delta(T*) ~ 1.4 k_alpha T* Pi^(-1/3) for Pi >> 1 (result of 6 Oct: ODE 4.3e-5
    against 1.4 x 3.0e-5 for the Table 2 value). The "law" column below is the
    t^(1/3) form and is not the relevant one on T* = 1.

    with eta cancelling in the balance k_alpha_bar alpha = mu c phi delta^2: the
    stress term stops the biofilm where the elastic energy mu c delta^2 reaches
    the growth energy k_alpha_bar = eta k_alpha. The script checks the law
    against the ODE and against the 3D seed surface of stress_to_growth_3d.py.

(B) A different growth law: Eq. 36 with a homeostatic pressure p_h,

        alpha_dot = k_alpha phi max(0, 1 - p / p_h),      p = -tr(sigma)/3 (compression > 0),

    and the phi equation as in Klempt 2024 Eq. 34 (no stress term). The
    biofilm stays; its volumetric growth stops where the pressure reaches p_h.
    Run on the 3D octant with the spherical seed (stress_to_growth_3d.Mesh),
    p_h relative to the largest seed pressure without feedback, p_0.

    python keio_wp2/homeostatic.py -> keio_wp2/results_homeostatic.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stress_to_growth_3d as S  # noqa: E402

K_ALPHA = S.K_ALPHA


def ode_A(mu_star, c, t_end=1.0):
    f = lambda t, y: [K_ALPHA * (1 + y[1]) - mu_star * c * y[0] * y[1] ** 2, K_ALPHA * y[0]]
    sol = solve_ivp(f, (0, t_end), [1.0, 0.0], method="Radau", rtol=1e-9, atol=1e-14, dense_output=True)
    t = np.linspace(0.01, t_end, 100)
    y = sol.sol(t)
    law = (3 * K_ALPHA ** 2 * t / (mu_star * c)) ** (1 / 3)
    return dict(t=t.tolist(), phi=y[0].tolist(), delta=y[1].tolist(), law=law.tolist())


def mechanics_p(M, phi, delta):
    """element mean pressure (compression > 0) [Pa]."""
    Ee = S.E0 * (phi ** 2 + S.FLOOR)
    mu = Ee / (2 * (1 + S.NU)); K = Ee / (3 * (1 - 2 * S.NU))
    vals = (mu[:, None, None] * M.Kdev + K[:, None, None] * M.Kvol).ravel()
    Kg = sp.csr_matrix((vals, (M.rows, M.cols)), shape=(M.ndof, M.ndof))
    f = np.zeros(M.ndof)
    np.add.at(f, M.dofs.ravel(), ((3 * K * delta)[:, None] * M.fvol).ravel())
    u = np.zeros(M.ndof)
    u[M.free] = spla.spsolve(Kg[M.free][:, M.free].tocsc(), f[M.free])
    e = (M.B0 @ u[M.dofs].T).T
    tr_e = e[:, 0] + e[:, 1] + e[:, 2] - 3 * delta
    return -K * tr_e


def run_B(n, p_h, dt=0.02, t_end=1.0):
    M = S.Mesh(n, "sphere")
    phi = M.inside.astype(float)
    alpha = np.ones_like(phi)
    hist = []
    for k in range(int(round(t_end / dt))):
        p = mechanics_p(M, phi, alpha - 1)
        g = np.ones_like(p) if p_h is None else np.clip(1 - p / p_h, 0.0, 1.0)
        phi = np.clip(phi + dt * (S.BETA * M.lap(phi) + K_ALPHA * alpha), 0, 1)
        alpha = alpha + dt * K_ALPHA * phi * g
        ins = M.inside
        hist.append(((k + 1) * dt, float((alpha - 1)[ins].mean()), float(phi[ins].mean()),
                     float(p[ins].mean()), float(p[ins].max()), float(g[ins].mean())))
    return dict(p_h=p_h, n=n, dt=dt, hist=hist)


def main():
    t0 = time.time()
    out = {"A": [], "B": []}
    d3 = json.loads((HERE / "results_3d_n16.json").read_text())
    for r in d3["runs"]:
        if r["seed"] != "sphere" or r["mu_star"] == 0:
            continue
        c = d3["runs"][0]["hist"][-1][6]                  # c at the surface without the term
        o = ode_A(r["mu_star"], c)
        o.update(mu_star=r["mu_star"], c=c)
        out["A"].append(o)
        print(f"A mu*={r['mu_star']:.3g}: ODE delta(1) {o['delta'][-1]:.3e}, law {o['law'][-1]:.3e}, "
              f"3D seed {r['hist'][-1][1]:.3e}", flush=True)
    base = run_B(16, None)
    out["B"].append(base)
    p0 = base["hist"][-1][4]
    print(f"B no feedback: alpha-1 {base['hist'][-1][1]:.3e}, p mean {base['hist'][-1][3]:.3e}, "
          f"max {p0:.3e} Pa [{time.time() - t0:.0f} s]", flush=True)
    for frac in (0.75, 0.5, 0.25):
        r = run_B(16, frac * p0)
        r["frac"] = frac
        out["B"].append(r)
        t, d, ph, pm, px, g = r["hist"][-1]
        print(f"B p_h = {frac} p0: alpha-1 {d:.3e}, phi {ph:.3f}, p mean {pm:.3e} max {px:.3e} "
              f"(p_h {frac * p0:.3e}), g {g:.2f} [{time.time() - t0:.0f} s]", flush=True)
    out["p0"] = p0
    (HERE / "results_homeostatic.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
