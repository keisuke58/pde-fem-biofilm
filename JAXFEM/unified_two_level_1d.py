"""unified_two_level_1d.py -- the two-level unified model (UNIFIED_MODEL_CONCEPT.ja.md, section 8),
1D, one residual for all fields, implicit Euler, full Newton.

Two levels, as the one-level prototype (unified_1d_monolithic.py) showed is needed:
  phi(x, t)            biofilm phase (Klempt et al. 2024): biofilm vs outer liquid
  s_i, s_0, psi_i      composition inside the biofilm (Klempt et al. 2026 point model):
                       species i and inner liquid s_0, sum_i s_i + s_0 = 1, living fractions psi_i
  absolute amounts     phi_i = phi s_i,  living biomass phi sum_i s_i psi_i
  alpha                local expansion, F_g = alpha I

Energy density and dissipation (w = phi + eps; eps keeps the composition equations regular in the
outer liquid, where they carry no weight):
  E = beta/2 |phi'|^2 + lam * w * Psi_pm(s, s_0, psi) + w * beta_s/2 sum_i |s_i'|^2
      - k_alpha (alpha - 1) * phi * sum_i s_i psi_i
  D = eta_phi/2 phidot^2 + lam * w * D_pm(sdot, s0dot, psidot) + 1/2 alphadot^2
  Psi_pm, D_pm: the point model of Klempt et al. 2026 (as in core_hamilton_1d.residual, eta_i = 1).
Constraint w * gamma (sum s + s_0 - 1).

Stationarity conditions, divided by w where they carry it:
  s_i, s_0, psi_i : the point model, with the growth term k_alpha (alpha - 1) (phi / (lam w)) in (A sbar)_i / c,
                    and the composition diffusion - beta_s (w s_i')' / w
  phi             : eta_phi phidot - beta phi'' + lam Psi_pm + beta_s/2 |s'|^2 - k_alpha (alpha - 1) sum s psi = 0
  alpha           : alphadot = k_alpha phi sum_i s_i psi_i      (scheme (c) of the element, derived)
(the state dependence of the dissipation weight is not a force in the extended Hamilton principle,
which varies the dissipation with respect to the rates only.)

Limits: lam -> 0 and beta_s = 0 give the element's scheme (c): the composition runs on its own at
every point (the free point model) and alpha grows with phi times the living biomass.
Only public sources (Klempt et al. 2024 and 2026) and this repository's code.
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
NV = 14   # phi, s1..5, s0, psi1..5, gamma, alpha
TOL = 1e-8  # residual tolerance; the composition rows in the outer liquid are divided by lam*w ~ 1e-7, so round-off there is ~1e-9
B = lambda x: 1.0 / (x ** 2 * (1.0 - x) ** 2)
dB = lambda x: KP1 * (2.0 - 4.0 * x) / ((x - 1.0) ** 3 * x ** 3)          # KP1 B'(x)


def unpack(U, N):
    X = U.reshape(N, NV)
    return X[:, 0], X[:, 1:6], X[:, 6], X[:, 7:12], X[:, 12], X[:, 13]


def psi_barrier(s, s0, psi):
    return KP1 * (jnp.sum(B(s), 1) + B(s0) + jnp.sum(B(psi), 1))


def psi_int(s, psi, P):
    sb = s * psi
    return 0.5 * P["ast"] * jnp.sum(P["b"] * psi ** 2, 1) - 0.5 * P["c"] * jnp.einsum("ki,ij,kj->k", sb, P["A"], sb)


def psi_pm(s, s0, psi, P):
    sb = s * psi
    return (KP1 * (jnp.sum(B(s), 1) + B(s0) + jnp.sum(B(psi), 1)) + 0.5 * P["ast"] * jnp.sum(P["b"] * psi ** 2, 1)
            - 0.5 * P["c"] * jnp.einsum("ki,ij,kj->k", sb, P["A"], sb))


def d_pm(s, s0, psi, sn, s0n, psin, dt):
    sd, s0d, pd = (s - sn) / dt, (s0 - s0n) / dt, (psi - psin) / dt
    return 0.5 * jnp.sum(sd ** 2 + (psi * sd + s * pd) ** 2, 1) + 0.5 * s0d ** 2


def grad_nodes(f, dx):
    """d/dx at cell midpoints"""
    return jnp.diff(f, axis=0) / dx


def weight(phi, P):
    """Weight of the composition level: h(phi) + eps.  "phi": h = phi (first version); "phi2":
    h = phi^2, as the mechanical term of Klempt et al. 2024 (phi^2 mu), so that dh/dphi = 0 in the
    outer liquid and the composition energy cannot create biofilm where there is none."""
    return (phi ** 2 if P.get("weight", "phi2") == "phi2" else phi) + P["eps"]


def functional(U, Un, P, F):
    """Discrete functional of one implicit Euler step: E_h(U) + dt * D_h + constraint, nodal
    quadrature hw.  The weights and coefficients of the dissipation and of the constraint are
    held fixed: the extended Hamilton principle varies D with respect to the rates only.  They
    are passed in F = (w, s, psi) at the new state, so that the gradient with respect to U treats
    them as constants, while the residual, as a function of U, still depends on them (its
    Jacobian, taken by jacfwd of residual(), includes that dependence).  The gradient of this
    functional with respect to U is the residual (see residual())."""
    N, dt, dx, eps = P["N"], P["dt"], P["dx"], P["eps"]
    Wf, Sf, Pf = F
    phi, s, s0, psi, gam, al = unpack(U, N)
    phin, sn, s0n, psin, _, aln = unpack(Un, N)
    hw = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    lam, ka, beta, bs = P["lam"], P["ka"], P["beta"], P["beta_s"]
    w = weight(phi, P)
    # the point model's barrier terms are a regularised indicator of the box constraints
    # 0 < s, s_0, psi < 1: their value is not stored energy (it is 517 at the calibration's start,
    # psi = 0.999, against 8 for psi = 0.5), so they act on the composition only and do no work
    # on the phase (weight held at w, like the constraint).  The physical energy Psi_int
    # (interactions, antibiotic) carries the weight w.
    # P["barrier"]: "frozen" = weight w held fixed (a regularised constraint, no force on phi);
    # "const" = constant weight lam (stored energy independent of phi, so the energy inequality
    # d/dt E = -2 D holds exactly, and still no force on phi)
    bmode = P.get("barrier", "stored")
    Wb = Wf if bmode == "frozen" else (1.0 if bmode == "const" else w)   # "stored": Klempt 2026 literal, weight w
    E = jnp.sum(hw * (lam * w * psi_int(s, psi, P) + lam * Wb * psi_barrier(s, s0, psi)
                      - ka * (al - 1.0) * phi * jnp.sum(s * psi, 1)))
    E += 0.5 * beta * jnp.sum(jnp.diff(phi) ** 2) / dx
    # phi <= 1 (a volume fraction): one-sided barrier, a regularised indicator like the point
    # model's barrier terms; KPHI small, so it acts only next to phi = 1
    E += jnp.sum(hw * P.get("kphi", 1e-5) / (1.0 - phi) ** 2)
    wm = 0.5 * (w[1:] + w[:-1])
    E += 0.5 * bs * jnp.sum(wm[:, None] * jnp.diff(s, axis=0) ** 2) / dx
    sd, pd, s0d = (s - sn) / dt, (psi - psin) / dt, (s0 - s0n) / dt
    S, Ps, W = Sf, Pf, Wf
    Dpm = 0.5 * jnp.sum(sd ** 2 + (Ps * sd + S * pd) ** 2, 1) + 0.5 * s0d ** 2
    D = jnp.sum(hw * (0.5 * P["eta_phi"] * ((phi - phin) / dt) ** 2 + lam * W * Dpm + 0.5 * ((al - aln) / dt) ** 2))
    C = jnp.sum(hw * lam * W * gam * (jnp.sum(s, 1) + s0 - 1.0))
    return E + dt * D + C


def energy(U, P):
    """Stored energy E_h of the state (the functional without dissipation and constraint)."""
    N, dx = P["N"], P["dx"]
    phi, s, s0, psi, gam, al = unpack(U, N)
    hw = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    w = weight(phi, P)
    bmode = P.get("barrier", "stored")
    Wb = 1.0 if bmode == "const" else w
    E = jnp.sum(hw * (P["lam"] * w * psi_int(s, psi, P) + P["lam"] * Wb * psi_barrier(s, s0, psi)
                      - P["ka"] * (al - 1.0) * phi * jnp.sum(s * psi, 1) + P.get("kphi", 1e-5) / (1.0 - phi) ** 2))
    E += 0.5 * P["beta"] * jnp.sum(jnp.diff(phi) ** 2) / dx
    wm = 0.5 * (w[1:] + w[:-1])
    return E + 0.5 * P["beta_s"] * jnp.sum(wm[:, None] * jnp.diff(s, axis=0) ** 2) / dx


def residual(U, Un, P):
    """Gradient of the discrete functional (JAX), rows scaled to the form of the point model:
    composition and constraint rows divided by lam * w * hw, phi and alpha rows by hw.  The
    equations are therefore the discrete Euler-Lagrange equations by construction."""
    N, dx, eps = P["N"], P["dx"], P["eps"]
    X = U.reshape(N, NV)
    w = weight(X[:, 0], P)
    G = jax.grad(functional)(U, Un, P, (w, X[:, 1:6], X[:, 7:12])).reshape(N, NV)
    hw = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    G = G.at[:, 1:13].divide((P["lam"] * w * hw)[:, None])
    G = G.at[:, 0].divide(hw).at[:, 13].divide(hw)
    return G.ravel()


def residual_hand(U, Un, P):
    N, dt, dx, eps = P["N"], P["dt"], P["dx"], P["eps"]
    phi, s, s0, psi, gam, al = unpack(U, N)
    phin, sn, s0n, psin, _, aln = unpack(Un, N)
    c, A, lam, ka, beta, bs = P["c"], P["A"], P["lam"], P["ka"], P["beta"], P["beta_s"]
    w = phi + eps
    sd, pd, s0d = (s - sn) / dt, (psi - psin) / dt, (s0 - s0n) / dt
    sb = s * psi
    Ia = sb @ A
    # growth energy -ka (al - 1) phi sum s psi ; its derivative w.r.t. s_i is -ka (al-1) phi psi_i
    g = ka * (al - 1.0) * phi / (lam * w)                  # per unit (lam w), enters like c (A sbar)_i
    # composition diffusion: energy w bs/2 |s'|^2 -> -(w_mid bs s')' ; divided by lam w below
    wm = 0.5 * (w[1:] + w[:-1])
    flux = bs * wm[:, None] * grad_nodes(s, dx)            # at midpoints
    div = jnp.zeros_like(s).at[:-1].add(flux).at[1:].add(-flux) / dx   # = (w s')' with no-flux ends (sign below)
    hw = jnp.full(N, dx).at[0].set(dx / 2).at[-1].set(dx / 2)
    div = div * dx / hw[:, None]
    q_s = (dB(s) + gam[:, None] + (1.0 + psi ** 2) * sd + s * psi * pd
           - c * psi * Ia - g[:, None] * psi - div / (lam * w)[:, None])
    q_s0 = gam + dB(s0) + s0d
    q_psi = (KP1 * ((-2.0) / ((psi - 1.0) ** 2 * psi ** 3) - 2.0 / ((psi - 1.0) ** 3 * psi ** 2))
             + P["b"] * P["ast"] * psi + s * psi * sd + s ** 2 * pd - c * s * Ia - g[:, None] * s)
    q_con = jnp.sum(s, 1) + s0 - 1.0
    # phi: eta_phi phidot - beta phi'' + lam (Psi_pm + D_pm-part) + bs/2 |s'|^2 - ka (al - 1) sum s psi = 0
    lapphi = jnp.zeros_like(phi)
    lapphi = lapphi.at[1:-1].set((phi[:-2] + phi[2:] - 2 * phi[1:-1]) / dx ** 2)
    lapphi = lapphi.at[0].set(2 * (phi[1] - phi[0]) / dx ** 2).at[-1].set(2 * (phi[-2] - phi[-1]) / dx ** 2)
    sg2 = jnp.zeros(N).at[:-1].add(0.5 * jnp.sum(grad_nodes(s, dx) ** 2, 1)).at[1:].add(0.5 * jnp.sum(grad_nodes(s, dx) ** 2, 1)) * 0.5
    q_phi = (P["eta_phi"] * (phi - phin) / dt - beta * lapphi
             + lam * psi_pm(s, s0, psi, P)
             + 0.5 * bs * sg2 - ka * (al - 1.0) * jnp.sum(sb, 1))
    q_al = (al - aln) / dt - ka * phi * jnp.sum(sb, 1)
    return jnp.concatenate([q_phi[:, None], q_s, q_s0[:, None], q_psi, q_con[:, None], q_al[:, None]], 1).ravel()


def make_newton(P):
    """Damped Newton for one implicit Euler step; dt is an argument so that the step can be
    subdivided without recompiling."""
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
            t = jnp.minimum(1.0, jnp.minimum(jnp.min(lim), jnp.min(limp)))
            U = U + t * d
            return U, k + 1, jnp.max(jnp.abs(R(U, Un, dt)))

        return jax.lax.while_loop(cond, body, (Un, 0, jnp.inf))
    return newton


def advance(newton, U, dt, depth=0, max_depth=12):
    """One step of length dt; if Newton does not converge, two half steps (recursively)."""
    V, it, r = newton(U, dt)
    if float(r) <= TOL:
        return V, int(it), float(r), 1
    if depth >= max_depth:
        raise RuntimeError(f"Newton failed at dt = {dt:.3e} (residual {float(r):.2e})")
    V, i1, r1, n1 = advance(newton, U, dt / 2, depth + 1, max_depth)
    V, i2, r2, n2 = advance(newton, V, dt / 2, depth + 1, max_depth)
    return V, max(i1, i2), max(r1, r2), n1 + n2


def run(phi0_x, g_comp, theta, *, lam, ka, beta, beta_s, eta_phi=1.0, eps=1e-6, ast=0.0, ab_on=None, weight="phi2", barrier="stored",
        T=0.2375, steps=2375, L=1.0, n_out=10):
    N = len(phi0_x)
    A, _ = theta_to_matrices(jnp.asarray(theta))
    P0 = {"N": N, "dt": T / steps, "dx": L / (N - 1), "c": C_STAR, "A": A, "b": jnp.ones(5), "ast": 0.0,
          "lam": lam, "ka": ka, "beta": beta, "beta_s": beta_s, "eta_phi": eta_phi, "eps": eps, "weight": weight, "barrier": barrier}
    P1 = dict(P0, ast=ast)
    n0, n1 = make_newton(P0), make_newton(P1)
    X = np.zeros((N, NV)); X[:, 0] = phi0_x
    s = np.clip(np.asarray(g_comp[0:5]), 1e-10, None); X[:, 1:6] = s; X[:, 6] = max(1 - s.sum(), 1e-10)
    X[:, 1:6] *= ((1 - X[:, 6]) / X[:, 1:6].sum(1))[:, None]
    X[:, 7:12] = np.clip(g_comp[6:11], 1e-10, 1 - 1e-10); X[:, 13] = 1.0
    U = jnp.asarray(X.ravel()); outs = [np.asarray(U)]; worst = 0.0; maxit = 0; nsteps_total = 0
    Es = [float(energy(U, P0))]
    save = set(int(round(steps * j / n_out)) for j in range(1, n_out + 1))
    for k in range(steps):
        on = ast > 0 and ab_on is not None and (k + 0.5) / steps >= ab_on
        U, it, r, nsub = advance(n1 if on else n0, U, T / steps)
        maxit = max(maxit, it); worst = max(worst, r); nsteps_total += nsub
        Es.append(float(energy(U, P1 if on else P0)))
        if k + 1 in save:
            outs.append(np.asarray(U))
    run.energy = np.array(Es)
    return np.array(outs), maxit, worst, nsteps_total


if __name__ == "__main__":
    import unified_1d as u1
    import unified_1d_monolithic as m1
    g0, theta = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1
    # check: lam -> small, beta_s = 0, uniform phi: the composition is the point model at every node
    outs, it, r, _ = run(np.full(3, 0.8), g, theta, lam=1e-6, ka=0.0, beta=0.0, beta_s=0.0, steps=200, T=0.02, n_out=1)
    ref = m1.point_model_converged(g, theta, 1e-4, 200)
    X = outs[-1].reshape(3, NV)
    print(f"check (uniform, lam = 1e-6, ka = 0): max |composition - point model| = {np.abs(X[0, 1:12] - ref[0:11]).max():.2e}, "
          f"phi {X[0, 0]:.6f}, Newton <= {it}, residual <= {r:.1e}")
    # seed column: phi = 0.9 at the substrate, outer liquid elsewhere
    x = np.linspace(0, 1, 41); phi_x = np.where(x <= 0.25, 0.9, 1e-4)
    for name, ast, on in (("no antibiotic", 0.0, None), ("antibiotic from T* = 0.5", 0.65, 0.5)):
        outs, it, r, ns = run(phi_x, g, theta, lam=1e-3, ka=2.0, beta=0.02, beta_s=0.02, ast=ast, ab_on=on, steps=119)
        X = outs[-1].reshape(41, NV)
        front = x[np.where(X[:, 0] > 0.45)[0].max()]
        print(f"{name:26s}: alpha(0) = {X[0, 13]:.4f}, phi at x = 0, 0.5, 1: {X[0, 0]:.3f} {X[20, 0]:.4f} {X[-1, 0]:.4f}, "
              f"front {front:.3f}, living share(0) = {(X[0, 1:6] * X[0, 7:12]).sum():.3f}, Newton <= {it}, residual <= {r:.1e}")
