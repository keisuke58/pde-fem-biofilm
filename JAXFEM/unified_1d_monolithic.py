"""unified_1d_monolithic.py -- the unified model of UNIFIED_MODEL_CONCEPT.ja.md solved as one
variational system: all fields of all nodes in one residual, implicit Euler, full Newton.

Unknowns per node: phi_1..5, phi_0, psi_1..5, gamma, alpha  (13).
Energy density (eta_i = 1, Kp1 = 1e-4, c = c* constant; Klempt et al. 2026 point model):
  E = sum_i Kp1 B(phi_i) + Kp1 B(phi_0) + sum_i Kp1 B(psi_i) + sum_i (b_i a*/2) psi_i^2
      - (c/2) phibar^T A phibar + sum_i (beta/2) |phi_i'|^2 - k_alpha (alpha - 1) sum_i phibar_i,
  B(x) = 1 / (x^2 (1 - x)^2),  phibar_i = phi_i psi_i.
Dissipation (Klempt et al. 2026 for phi, psi, phi_0; Klempt et al. 2024 for alpha):
  D = sum_i [1/2 phidot_i^2 + 1/2 (d phibar_i/dt)^2] + 1/2 phi0dot^2 + 1/2 alphadot^2.
Constraint gamma (sum_i phi_i + phi_0 - 1).  The stationarity conditions are
  dE/dphi_i - beta phi_i'' + dD/dphidot_i + gamma = 0, dE/dphi_0 + phi0dot + gamma = 0,
  dE/dpsi_i + dD/dpsidot_i = 0, alphadot - k_alpha sum_i phibar_i = 0, sum phi + phi_0 = 1,
with n . grad phi_i = 0 at both ends.  At a node with uniform neighbours and k_alpha = 0 they are
exactly the point model of core_hamilton_1d.residual (gate off).  Left out on purpose because they
are not variational: Klempt 2024's front term and the nutrient field (c = c*).

Because the evolution is a gradient flow of E with the quadratic D, E must not increase along the
solution (no external power).  The check below reports the largest increase of the discrete E.
Only public sources (Klempt et al. 2024, 2026) and this repository's code are used.

Run:  python JAXFEM/unified_1d_monolithic.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "JAXFEM")); sys.path.insert(0, str(ROOT / "ansys_usermat/coupling"))
from core_hamilton_1d import theta_to_matrices  # noqa: E402
from ecology_constants import C_STAR  # noqa: E402

KP1 = 1e-4
NV = 13


def unpack(U, N):
    X = U.reshape(N, NV)
    return X[:, 0:5], X[:, 5], X[:, 6:11], X[:, 11], X[:, 12]


def lap(phi, dx):
    """second difference with mirrored ends (n . grad phi = 0), the variation of the
    discrete gradient energy sum_k beta/2 ((phi_{k+1} - phi_k)/dx)^2 dx, divided by the nodal weight"""
    l = jnp.zeros_like(phi)
    l = l.at[1:-1].set((phi[:-2] + phi[2:] - 2 * phi[1:-1]) / dx ** 2)
    l = l.at[0].set(2 * (phi[1] - phi[0]) / dx ** 2).at[-1].set(2 * (phi[-2] - phi[-1]) / dx ** 2)
    return l


def residual(U, Un, P):
    N, dt, dx = P["N"], P["dt"], P["dx"]
    phi, phi0, psi, gam, al = unpack(U, N)
    phin, phi0n, psin, _, aln = unpack(Un, N)
    c, A, b, ast, ka, beta = P["c"], P["A"], P["b"], P["ast"], P["ka"], P["beta"]
    phid = (phi - phin) / dt; psid = (psi - psin) / dt; phi0d = (phi0 - phi0n) / dt
    Ia = (phi * psi) @ A + (ka * (al - 1.0) / c)[:, None]
    q_phi = (KP1 * (2.0 - 4.0 * phi) / ((phi - 1.0) ** 3 * phi ** 3)
             + gam[:, None] + (1.0 + psi ** 2) * phid + phi * psi * psid
             - c * psi * Ia - beta * lap(phi, dx))
    q_phi0 = gam + KP1 * (2.0 - 4.0 * phi0) / ((phi0 - 1.0) ** 3 * phi0 ** 3) + phi0d
    q_psi = ((-2.0 * KP1) / ((psi - 1.0) ** 2 * psi ** 3) - (2.0 * KP1) / ((psi - 1.0) ** 3 * psi ** 2)
             + b * ast * psi + phi * psi * phid + phi ** 2 * psid - c * phi * Ia)
    q_con = jnp.sum(phi, axis=1) + phi0 - 1.0
    q_al = (al - aln) / dt - ka * jnp.sum(phi * psi, axis=1)
    return jnp.concatenate([q_phi, q_phi0[:, None], q_psi, q_con[:, None], q_al[:, None]], axis=1).ravel()


def energy(U, P):
    N, dx = P["N"], P["dx"]
    phi, phi0, psi, _, al = unpack(U, N)
    B = lambda x: 1.0 / (x ** 2 * (1.0 - x) ** 2)
    pb = phi * psi
    e = (KP1 * (jnp.sum(B(phi), 1) + B(phi0) + jnp.sum(B(psi), 1))
         + 0.5 * P["ast"] * jnp.sum(P["b"] * psi ** 2, 1)
         - 0.5 * P["c"] * jnp.einsum("ki,ij,kj->k", pb, P["A"], pb)
         - P["ka"] * (al - 1.0) * jnp.sum(pb, 1))
    w = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    grad = 0.5 * P["beta"] * jnp.sum(jnp.diff(phi, axis=0) ** 2, axis=1).sum() / dx
    return jnp.sum(w * e) + grad


def make_step(P):
    R = lambda U, Un: residual(U, Un, P)
    J = jax.jacfwd(R, argnums=0)

    @jax.jit
    def newton(Un):
        def cond(s):
            U, k, r = s
            return (r > 1e-11) & (k < 30)

        def body(s):
            U, k, _ = s
            F = R(U, Un)
            d = jnp.linalg.solve(J(U, Un), -F)
            # damping: keep phi, phi_0, psi inside (0, 1)
            X = U.reshape(-1, NV); D = d.reshape(-1, NV)
            idx = jnp.r_[0:11]
            x, dd = X[:, idx], D[:, idx]
            lim = jnp.where(dd < 0, -0.9 * x / jnp.minimum(dd, -1e-300), jnp.where(dd > 0, 0.9 * (1 - x) / jnp.maximum(dd, 1e-300), jnp.inf))
            t = jnp.minimum(1.0, jnp.min(lim))
            U = U + t * d
            return U, k + 1, jnp.max(jnp.abs(R(U, Un)))

        U, k, r = jax.lax.while_loop(cond, body, (Un, 0, jnp.inf))
        return U, k, r
    return newton


def build(theta, N, *, ka, beta, ast=0.0, b=None, T=0.2375, steps=2375, L=1.0):
    A, _ = theta_to_matrices(jnp.asarray(theta))
    P = {"N": N, "dt": T / steps, "dx": L / (N - 1), "c": C_STAR, "A": A,
         "b": jnp.ones(5) if b is None else jnp.asarray(b, float), "ast": ast, "ka": ka, "beta": beta}
    return P


def initial(g_nodes):
    """g_nodes (N, 12) point-model states -> U with alpha = 1; gamma made consistent by one solve
    is not needed: gamma is algebraic and set by the first Newton step."""
    N = g_nodes.shape[0]
    X = np.zeros((N, NV)); X[:, 0:12] = g_nodes; X[:, 12] = 1.0
    X[:, 0:11] = np.clip(X[:, 0:11], 1e-10, 1 - 1e-10)      # as core_hamilton_1d.clip_state
    X[:, 5] = np.clip(1.0 - X[:, 0:5].sum(1), 1e-10, None)
    X[:, 0:5] *= ((1.0 - X[:, 5]) / X[:, 0:5].sum(1))[:, None]
    return jnp.asarray(X.ravel())


def run(g_nodes, theta, *, ka, beta, ast=0.0, ab_on=None, b=None, T=0.2375, steps=2375, L=1.0, n_out=10):
    N = g_nodes.shape[0]
    P0 = build(theta, N, ka=ka, beta=beta, ast=0.0, b=b, T=T, steps=steps, L=L)
    P1 = dict(P0, ast=ast)
    step0, step1 = make_step(P0), make_step(P1)
    U = initial(g_nodes)
    E = [float(energy(U, P0))]; outs = [np.asarray(U)]; worst = 0.0; maxit = 0
    save = set(int(round(steps * j / n_out)) for j in range(1, n_out + 1))
    for k in range(steps):
        on = ast > 0 and ab_on is not None and (k + 0.5) / steps >= ab_on
        P = P1 if on else P0
        Un = U
        U, it, r = (step1 if on else step0)(Un)
        maxit = max(maxit, int(it)); worst = max(worst, float(r))
        e = float(energy(U, P)); E.append(e)
        if k + 1 in save:
            outs.append(np.asarray(U))
    return np.array(outs), np.array(E), maxit, worst


def point_model_converged(g, theta, dt, steps):
    """Independent reference: core_hamilton_1d.residual (gate off) solved to |Q| < 1e-11 every
    step with damped Newton.  (ecology_jax's fixed six iterations leave the first ~50 steps
    from the calibration's start state unconverged; see UNIFIED_MODEL_CONCEPT.ja.md.)"""
    import core_hamilton_1d as ch
    import ecology_jax as ej
    A, b = theta_to_matrices(jnp.asarray(theta))
    p = jax.tree_util.tree_map(jnp.asarray, ej.default_hparams(dt)); p["A"] = A; p["b_diag"] = b
    F = jax.jit(lambda gn, gp: ch.residual(gn, gp, p)); J = jax.jit(jax.jacfwd(lambda gn, gp: ch.residual(gn, gp, p)))
    g = jnp.asarray(g)
    for _ in range(steps):
        gp, x = g, g
        for _ in range(40):
            r = F(x, gp)
            if float(jnp.max(jnp.abs(r))) < 1e-11:
                break
            d = jnp.linalg.solve(J(x, gp), -r)
            v, dv = x[0:11], d[0:11]
            lim = jnp.where(dv < 0, -0.9 * v / jnp.minimum(dv, -1e-300), jnp.where(dv > 0, 0.9 * (1 - v) / jnp.maximum(dv, 1e-300), jnp.inf))
            x = x + min(1.0, float(jnp.min(lim))) * d
        g = x
    return np.asarray(g)


if __name__ == "__main__":
    import unified_1d as u1
    g0, theta = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1        # interior start (phi_0 = 0.1)
    # check 1: three identical nodes, ka = 0, beta > 0: every node is the point model (converged implicit Euler)
    outs, E, it, r = run(np.repeat(g[None], 3, 0), theta, ka=0.0, beta=0.02, steps=200, T=0.02, n_out=1)
    gref = point_model_converged(g, theta, 1e-4, 200)
    X = outs[-1].reshape(3, NV)
    print(f"check 1 (uniform, ka = 0): max |node - converged point model| = {np.abs(X[0, 0:12] - gref).max():.2e}, "
          f"Newton <= {it} iterations, residual <= {r:.1e}, largest energy increase {np.diff(E).max():.2e}")
    # check 2: seed column; the void starts at the barrier floor phi = 0.05 (the barrier energy keeps every phi_i > 0)
    x, G = u1.seed_column(g0, N=41, phi_out=0.05)
    for name, ast, on in (("no antibiotic", 0.0, None), ("antibiotic from T* = 0.5", 0.65, 0.5)):
        outs, E, it, r = run(G, theta, ka=2.0, beta=0.02, ast=ast, ab_on=on)
        X = outs[-1].reshape(41, NV)
        tot = X[:, 0:5].sum(1)
        print(f"{name:26s}: alpha(0) = {X[0, 12]:.4f}, sum phi at x = 0, 0.5, 1: {tot[0]:.3f} {tot[20]:.3f} {tot[-1]:.3f}, "
              f"living biomass(0) = {(X[0, 0:5] * X[0, 6:11]).sum():.3f}, Newton <= {it}, residual <= {r:.1e}, "
              f"largest energy increase {np.diff(E).max():.2e}")
