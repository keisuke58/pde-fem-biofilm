"""unified_two_level_fuel_1d.py -- the two-level model with growth as a nutrient-consuming reaction
(Onsager's variational principle; UNIFIED_MODEL_CONCEPT.ja.md section 9).

Starts from unified_two_level_1d.py and changes only the growth:
  * the growth energy -k_alpha (alpha - 1) phi sum s psi is removed (k_alpha = 0 there), so the
    growth no longer acts on psi (the side effect of section 8);
  * a nutrient field c with the free energy of an ideal solute, f(c) = c (ln c - 1), mu_c = ln c;
  * growth is a reaction that turns Y units of nutrient into one unit of growth xi = alpha - 1,
    with the chemical potential mu_b of the produced biomass.  Onsager: minimise the Rayleighian
        dE/dt + xidot^2 / (2 L) + |j|^2 / (2 D c)
    under the balance  cdot + Y xidot + div j = 0.  This gives
        xidot = L (Y ln c - mu_b),     j = -D grad c,     cdot = D lap c - Y xidot.
  * the living biomass is the catalyst: it sits in the kinetic coefficient, L = L0 phi sum_i s_i psi_i,
    not in the energy.  Lysis or an antibiotic lowers L and so the growth, with no force on psi.
Energy: the reaction lowers E by (Y mu_c - mu_b) xidot >= 0 per unit time wherever it runs forward;
with a closed boundary the total energy (two-level part + nutrient) must not increase.  An open
top boundary (c = 1, a reservoir) supplies nutrient, i.e. external work.
Only public sources (Klempt et al. 2024, 2026; Onsager's principle) and this repository's code.
"""
from __future__ import annotations

import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "JAXFEM", ROOT / "ansys_usermat/coupling"):
    sys.path.insert(0, str(p))
import unified_two_level_1d as t2  # noqa: E402
from core_hamilton_1d import theta_to_matrices  # noqa: E402
from ecology_constants import C_STAR  # noqa: E402

NV = t2.NV + 1          # ... , alpha (13), c (14)
TOL = 1e-8


def residual(U, Un, P):
    N, dt, dx = P["N"], P["dt"], P["dx"]
    X, Xn = U.reshape(N, NV), Un.reshape(N, NV)
    R14 = t2.residual(X[:, :t2.NV].ravel(), Xn[:, :t2.NV].ravel(), dict(P, ka=0.0)).reshape(N, t2.NV)
    phi, s, psi, al, c = X[:, 0], X[:, 1:6], X[:, 7:12], X[:, 13], X[:, 14]
    aln, cn = Xn[:, 13], Xn[:, 14]
    L = P["L0"] * phi * jnp.sum(s * psi, 1)
    xid = (al - aln) / dt
    r_al = xid - L * (P["Y"] * jnp.log(c) - P["mu_b"])
    lap = jnp.zeros(N)
    lap = lap.at[1:-1].set((c[:-2] + c[2:] - 2 * c[1:-1]) / dx ** 2)
    lap = lap.at[0].set(2 * (c[1] - c[0]) / dx ** 2).at[-1].set(2 * (c[-2] - c[-1]) / dx ** 2)
    r_c = (c - cn) / dt - P["Dc"] * lap + P["Y"] * xid
    if P["open_top"]:
        r_c = r_c.at[-1].set(c[-1] - 1.0)          # reservoir at the top
    R = jnp.concatenate([R14[:, :13], r_al[:, None], r_c[:, None]], 1)
    return R.ravel()


def total_energy(U, P):
    """two-level energy (barrier terms included as stored energy, growth energy off) + nutrient + mu_b xi"""
    N, dx = P["N"], P["dx"]
    X = U.reshape(N, NV)
    phi, s, s0, psi, al, c = X[:, 0], X[:, 1:6], X[:, 6], X[:, 7:12], X[:, 13], X[:, 14]
    hw = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    w = t2.weight(phi, P)
    e = P["lam"] * w * t2.psi_pm(s, s0, psi, P) + P.get("kphi", 1e-5) / (1 - phi) ** 2 + c * (jnp.log(c) - 1) + P["mu_b"] * (al - 1)
    E = jnp.sum(hw * e) + 0.5 * P["beta"] * jnp.sum(jnp.diff(phi) ** 2) / dx
    wm = 0.5 * (w[1:] + w[:-1])
    return E + 0.5 * P["beta_s"] * jnp.sum(wm[:, None] * jnp.diff(s, axis=0) ** 2) / dx


def make_newton(P):
    def R(U, Un, dt):
        return residual(U, Un, dict(P, dt=dt))
    J = jax.jacfwd(R, argnums=0)
    box = jnp.r_[1:12]

    @jax.jit
    def newton(Un, dt):
        def cond(st):
            return (st[2] > TOL) & (st[1] < 40)

        def body(st):
            U, k, _ = st
            d = jnp.linalg.solve(J(U, Un, dt), -R(U, Un, dt))
            X, D = U.reshape(-1, NV), d.reshape(-1, NV)
            x, dd = X[:, box], D[:, box]
            lim = jnp.where(dd < 0, -0.9 * x / jnp.minimum(dd, -1e-300), jnp.where(dd > 0, 0.9 * (1 - x) / jnp.maximum(dd, 1e-300), jnp.inf))
            ph, dph = X[:, 0], D[:, 0]
            limp = jnp.where(dph < 0, -0.9 * ph / jnp.minimum(dph, -1e-300), 0.9 * (1 - ph) / jnp.maximum(dph, 1e-300))
            cc, dc = X[:, 14], D[:, 14]
            limc = jnp.where(dc < 0, -0.9 * cc / jnp.minimum(dc, -1e-300), jnp.inf)
            t = jnp.minimum(1.0, jnp.minimum(jnp.min(lim), jnp.minimum(jnp.min(limp), jnp.min(limc))))
            U = U + t * d
            return U, k + 1, jnp.max(jnp.abs(R(U, Un, dt)))

        return jax.lax.while_loop(cond, body, (Un, 0, jnp.inf))
    return newton


def run(phi0_x, g_comp, theta, *, lam=1.0, beta=0.02, beta_s=0.02, L0=0.7, Y=1.0, c_stop=0.1, Dc=1.0,
        open_top=True, ast=0.0, ab_on=None, T=0.2375, steps=60, L=1.0, n_out=10, c0=1.0):
    N = len(phi0_x)
    A, _ = theta_to_matrices(jnp.asarray(theta))
    P0 = {"N": N, "dt": T / steps, "dx": L / (N - 1), "c": C_STAR, "A": A, "b": jnp.ones(5), "ast": 0.0,
          "lam": lam, "ka": 0.0, "beta": beta, "beta_s": beta_s, "eta_phi": 1.0, "eps": 1e-6, "weight": "phi2",
          "L0": L0, "Y": Y, "mu_b": Y * np.log(c_stop), "Dc": Dc, "open_top": open_top}
    P1 = dict(P0, ast=ast)
    n0, n1 = make_newton(P0), make_newton(P1)
    X = np.zeros((N, NV)); X[:, 0] = phi0_x
    s = np.clip(np.asarray(g_comp[0:5]), 1e-10, None); X[:, 1:6] = s; X[:, 6] = max(1 - s.sum(), 1e-10)
    X[:, 1:6] *= ((1 - X[:, 6]) / X[:, 1:6].sum(1))[:, None]
    X[:, 7:12] = np.clip(g_comp[6:11], 1e-10, 1 - 1e-10); X[:, 13] = 1.0; X[:, 14] = c0
    U = jnp.asarray(X.ravel()); outs = [np.asarray(U)]; E = [float(total_energy(U, P0))]
    save = set(int(round(steps * j / n_out)) for j in range(1, n_out + 1))

    def adv(newton, U, dt, depth=0):
        V, it, r = newton(U, dt)
        if float(r) <= TOL:
            return V
        if depth >= 12:
            raise RuntimeError(f"Newton failed at dt = {dt:.3e} (residual {float(r):.2e})")
        return adv(newton, adv(newton, U, dt / 2, depth + 1), dt / 2, depth + 1)

    for k in range(steps):
        on = ast > 0 and ab_on is not None and (k + 0.5) / steps >= ab_on
        U = adv(n1 if on else n0, U, T / steps)
        E.append(float(total_energy(U, P1 if on else P0)))
        if k + 1 in save:
            outs.append(np.asarray(U))
    return np.array(outs), np.array(E)


if __name__ == "__main__":
    import unified_1d as u1
    g0, th = u1.ch_case(); g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1
    x = np.linspace(0, 1, 41); phi_x = np.where(x <= 0.25, 0.9, 1e-4)
    # closed system: the total energy must not increase
    outs, E = run(phi_x, g, th, open_top=False, steps=30)
    print(f"closed boundary: largest energy increase {np.diff(E).max():.2e} (E {E[0]:.4f} -> {E[-1]:.4f})")
    for name, ast, on in (("no antibiotic", 0.0, None), ("antibiotic from T* = 0.5", 0.65, 0.5)):
        outs, E = run(phi_x, g, th, ast=ast, ab_on=on)
        X = outs[-1].reshape(41, NV)
        print(f"{name:26s}: alpha-1 at x = 0, 0.2: {X[0, 13] - 1:.4f} {X[8, 13] - 1:.4f}, c at x = 0, 0.25, 1: "
              f"{X[0, 14]:.3f} {X[10, 14]:.3f} {X[-1, 14]:.3f}, phi(0) {X[0, 0]:.3f}, living(0) {(X[0, 1:6] * X[0, 7:12]).sum():.3f}")
