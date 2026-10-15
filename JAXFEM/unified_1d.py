"""unified_1d.py -- 1D prototype of the unified model of UNIFIED_MODEL_CONCEPT.ja.md.

Every species is a field phi_i(x, t), psi_i(x, t) with one void phi_0 = 1 - sum phi_i.
Energy = point model of Klempt et al. 2026 (as in core_hamilton_1d.residual)
       + sum_i beta/2 |grad phi_i|^2              (Klempt et al. 2024 Eq. 20, per species)
       - k_alpha * alpha * sum_i phi_i psi_i       (Klempt 2024's -k_alpha alpha phi, with the living biomass)
The growth energy has the same structure as the interaction energy -(c/2) phibar^T A phibar:
its derivatives are -k_alpha alpha psi_i (w.r.t. phi_i) and -k_alpha alpha phi_i (w.r.t. psi_i),
so in the point-model residual it enters as (A phibar)_i + k_alpha alpha / c.  The alpha
equation is alpha_dot = k_alpha sum_i phi_i psi_i (eta_alpha = 1).

Solved by operator splitting per step: point-model Newton at every node with the growth term,
explicit diffusion of every phi_i (no-flux ends), explicit alpha update.  Splitting is a
numerical choice of the prototype, not part of the model.  Only public sources are used
(Klempt et al. 2024 and 2026, this repository's point-model residual); nothing from
F. Klempt's unpublished two-species work.

Checks (python JAXFEM/unified_1d.py):
  1. uniform field, k_alpha = 0, beta = 0  -> every node equals the point model alone (ecology_jax);
  2. seed at the substrate, CH MAP: alpha and the biofilm extent with and without an
     antibiotic switched on at half time (assumed b_i = 1).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "JAXFEM")); sys.path.insert(0, str(ROOT / "ansys_usermat/coupling"))
from core_hamilton_1d import theta_to_matrices, clip_state  # noqa: E402
from ecology_constants import C_STAR  # noqa: E402


def residual(g_new, g_prev, kag, p):
    """Point-model residual of core_hamilton_1d.residual (gate off) with the growth energy:
    kag = k_alpha * alpha at this node.  Same arithmetic otherwise."""
    dt, Kp1, Eta, EtaPhi, c, ast, A, b = p["dt"], 1e-4, jnp.ones(5), jnp.ones(5), p["c"], p["ast"], p["A"], p["b"]
    phi, phi0, psi, gam = g_new[0:5], g_new[5], g_new[6:11], g_new[11]
    phid = (phi - g_prev[0:5]) / dt; phi0d = (phi0 - g_prev[5]) / dt; psid = (psi - g_prev[6:11]) / dt
    Ia = A @ (phi * psi) + kag / c
    t1 = Kp1 * (2.0 - 4.0 * phi) / ((phi - 1.0) ** 3 * phi ** 3)
    t2 = (1.0 / Eta) * (gam + (EtaPhi + Eta * psi ** 2) * phid + Eta * phi * psi * psid)
    t3 = (c / Eta) * psi * Ia
    q_phi = t1 + t2 - t3
    q_phi0 = gam + Kp1 * (2.0 - 4.0 * phi0) / ((phi0 - 1.0) ** 3 * phi0 ** 3) + phi0d
    s1 = (-2.0 * Kp1) / ((psi - 1.0) ** 2 * psi ** 3) - (2.0 * Kp1) / ((psi - 1.0) ** 3 * psi ** 2)
    s2 = (b * ast / Eta) * psi
    s3 = phi * psi * phid + phi ** 2 * psid
    s4 = (c / Eta) * phi * Ia
    q_psi = s1 + s2 + s3 - s4
    return jnp.concatenate([q_phi, jnp.array([q_phi0]), q_psi, jnp.array([jnp.sum(phi) + phi0 - 1.0])])


def newton(g_prev, kag, p):
    mask = jnp.ones(5, dtype=jnp.int64)

    def body(g, _):
        g = clip_state(g, mask)
        F = lambda gg: residual(gg, g_prev, kag, p)
        return clip_state(g + jnp.linalg.solve(jax.jacfwd(F)(g), -F(g)), mask), None

    g, _ = jax.lax.scan(body, clip_state(g_prev, mask), jnp.arange(6))
    return g


newton_nodes = jax.jit(jax.vmap(newton, in_axes=(0, 0, None)))


def diffuse(G, beta, dt, dx):
    phi = G[:, 0:5]
    lap = jnp.zeros_like(phi)
    lap = lap.at[1:-1].set((phi[:-2] + phi[2:] - 2 * phi[1:-1]) / dx ** 2)
    lap = lap.at[0].set(2 * (phi[1] - phi[0]) / dx ** 2).at[-1].set(2 * (phi[-2] - phi[-1]) / dx ** 2)
    phi = jnp.clip(phi + dt * beta * lap, 1e-10, 1.0)
    tot = jnp.minimum(phi.sum(1), 1 - 1e-10)
    phi = phi * (tot / phi.sum(1))[:, None]
    return G.at[:, 0:5].set(phi).at[:, 5].set(1 - tot)


def ch_case():
    pm = json.load(open(ROOT / "ansys_usermat/apdl/results/2026-10-cal5/w8_CH_cal_k0882_pm.json"))
    C = {k: i for i, k in enumerate(pm["columns"])}
    rows = [r for r in pm["seed_trace"] if r[C["ip"]] == pm["seed_trace"][0][C["ip"]]]
    r0 = min(rows, key=lambda r: int(r[C["isubst"]]))
    g0 = np.array([float(r0[C[f"gold_{i}"]]) for i in range(1, 13)])
    return g0, np.array([float(x) for x in pm["theta"]])


def simulate(g_nodes, theta, *, k_alpha, beta, T=0.2375, steps=2375, L=1.0, ab=0.0, ab_on=0.5, b=None, shift=0.0):
    """g_nodes (N, 12) initial states; returns times, states (n_out, N, 12), alpha (n_out, N).
    shift = 0: growth energy -k_alpha alpha sum phibar (Klempt 2024's form, a source k_alpha alpha psi_i
    in every phi_i equation, also in the void); shift = 1: -k_alpha (alpha - 1) sum phibar (no source
    until growth has happened).  The alpha equation is the same in both."""
    A, b0 = theta_to_matrices(jnp.asarray(theta))
    b = jnp.ones(5) if b is None else jnp.asarray(b, float)
    N = g_nodes.shape[0]; dx = L / (N - 1); dt = T / steps
    p_off = {"dt": dt, "c": C_STAR, "ast": 0.0, "A": A, "b": b}
    p_on = dict(p_off, ast=ab)
    G = jnp.asarray(g_nodes); alpha = jnp.ones(N)
    out_t, out_G, out_a = [0.0], [np.asarray(G)], [np.asarray(alpha)]
    save = set(int(round(steps * j / 10)) for j in range(1, 11))
    for k in range(steps):
        p = p_on if (ab > 0 and (k + 0.5) / steps >= ab_on) else p_off
        G = newton_nodes(G, k_alpha * (alpha - shift), p)
        if beta > 0:
            G = diffuse(G, beta, dt, dx)
        alpha = alpha + dt * k_alpha * jnp.sum(G[:, 0:5] * G[:, 6:11], axis=1)
        if (k + 1) in save:
            out_t.append((k + 1) / steps); out_G.append(np.asarray(G)); out_a.append(np.asarray(alpha))
    return np.array(out_t), np.array(out_G), np.array(out_a)


def seed_column(g0, N=41, seed=0.25, phi_seed=0.9, phi_out=1e-8):
    x = np.linspace(0, 1, N)
    G = np.zeros((N, 12))
    shares = g0[0:5] / g0[0:5].sum()
    for j, xx in enumerate(x):
        tot = phi_seed if xx <= seed else phi_out
        G[j, 0:5] = shares * tot; G[j, 5] = 1 - tot; G[j, 6:11] = g0[6:11]
    return x, G


if __name__ == "__main__":
    import ecology_jax as ej
    g0, theta = ch_case()
    # check 1: uniform field, no growth coupling, no diffusion -> point model alone
    G = np.repeat(g0[None, :], 3, axis=0)
    t, Gs, _ = simulate(G, theta, k_alpha=0.0, beta=0.0, steps=2375)
    g_ref, _ = ej.ecology_substeps(g0, theta, 0.2375, 2375)
    print("check 1, max |node - point model| at T* = 1:", float(np.abs(Gs[-1, 0] - np.asarray(g_ref)).max()))
    # check 2: seed column, growth from the living biomass, antibiotic from half time
    x, G = seed_column(g0)
    res = {}
    for name, ab, sh in (("none, alpha", 0.0, 0.0), ("none, alpha-1", 0.0, 1.0), ("antibiotic, alpha-1", 0.65, 1.0)):
        t, Gs, a = simulate(G, theta, k_alpha=2.0, beta=0.02, ab=ab, shift=sh)
        tot = Gs[:, :, 0:5].sum(2)
        front = [x[np.where(tt > 0.45)[0].max()] if (tt > 0.45).any() else 0 for tt in tot]
        res[name] = (t, a[:, 0], front, (Gs[-1, :, 0:5] * Gs[-1, :, 6:11]).sum(1)[0])
        print(f"{name:20s}: alpha at the substrate {a[-1, 0]:.4f}, front (sum phi > 0.45) {front[0]:.3f} -> {front[-1]:.3f}, "
              f"living biomass at the substrate {res[name][3]:.3f}")
