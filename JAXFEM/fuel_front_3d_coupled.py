"""
Time (corrected 15 Oct 2026): the field and the point model share one time, T = 0.2375 = the
calibration's Day 1 to Day 21, as in the 1D models (s = 1).  The first runs used s = 0.2375, which
advanced the point model only to about Day 5.
fuel_front_3d_coupled.py -- step 2: the composition acts on the phase in 3D (lam > 0).

Two-level model of p3/unified_model.tex in 3D on a regular grid, beta_s = 0:
  composition  : the calibrated five-species point model at every node of the biofilm (with
                 beta_s = 0 its equations are the point model, independent of lam); nodes that
                 the phase reaches later start from the calibration's Day-1 state;
  phase        : eta_phi phidot = beta lap phi - 2 K_phi/(1 - phi)^3 - lam h'(phi) Psi_pm(s, psi),
                 h = phi^2, Psi_pm the point model's stored energy (barrier terms included);
  growth       : xidot = L0 phi sum_i s_i psi_i (Y ln c - mu_b),  cdot = D lap c - Y xidot;
  mechanics    : once at T* = 1 (fuel_front_3d.elasticity).
Splitting per step: point model (implicit, converged), phase (explicit), nutrient and growth
(explicit).  The splitting is a numerical choice; the 1D monolithic solver
(unified_two_level_1d.py) is the reference for the equations.
Cases: lam = 0 (the decoupled limit, as fuel_front_3d_species.py) and lam = 1, each without and
with an antibiotic from T* = 0.5.  Output: assets/fig_fuel_front_3d_coupled.pdf.
"""
import sys
from pathlib import Path

import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "JAXFEM", ROOT / "ansys_usermat/coupling", ROOT / "ansys_usermat"):
    sys.path.insert(0, str(p))
import core_hamilton_1d as ch  # noqa: E402
import fuel_front_3d as f3  # noqa: E402
import fuel_front_3d_species as f3s  # noqa: E402
import unified_1d as u1  # noqa: E402
from ecology_constants import C_STAR  # noqa: E402

KP1 = 1e-4


def psi_pm(G, A, ast, b, c=C_STAR):
    s, s0, psi = G[:, 0:5], G[:, 5], G[:, 6:11]
    B = lambda x: 1.0 / (x ** 2 * (1 - x) ** 2)
    sb = s * psi
    return (KP1 * (B(s).sum(1) + B(s0) + B(psi).sum(1)) + 0.5 * ast * (b * psi ** 2).sum(1)
            - 0.5 * c * np.einsum("ki,ij,kj->k", sb, A, sb))


def run(n=16, lam=1.0, beta=0.02, eta_phi=1.0, kphi=1e-5, D=0.2, L0=2.0, Y=1.0, c_stop=0.1,
        T=0.2375, s=1.0, ast=0.0, ab_on=0.5, b=None, phi_min=0.05, psi_ref=0.0):
    g0, th = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1; g[6:11] = np.clip(g[6:11], 1e-10, 1 - 1e-10)
    A = np.asarray(ch.theta_to_matrices(jnp.asarray(th))[0])
    b = np.ones(5) if b is None else np.asarray(b, float)
    h = 1.0 / n
    x = np.linspace(0, 1, n + 1)
    X, Yg, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt((X - 0.5) ** 2 + (Yg - 0.5) ** 2 + Z ** 2)
    phi = 0.5 * (1 - np.tanh((r - 0.25) / 0.03)) * 0.9 + 1e-4
    G = np.tile(g, (phi.size, 1))                       # composition at every node (used where phi > phi_min)
    started = (phi > phi_min).ravel()
    dt = 0.15 * h * h / max(D, beta / eta_phi); steps = int(np.ceil(T / dt)); dt = T / steps
    step0 = f3s.make_pm_step(th, s * dt, 0.0, b); step1 = f3s.make_pm_step(th, s * dt, ast, b)
    c = np.full_like(phi, c_stop); c[-1] = 1.0; xi = np.zeros_like(phi); mu_b = Y * np.log(c_stop)

    def lap(f):
        fp = np.pad(f, 1, mode="edge")
        return (fp[2:, 1:-1, 1:-1] + fp[:-2, 1:-1, 1:-1] + fp[1:-1, 2:, 1:-1] + fp[1:-1, :-2, 1:-1]
                + fp[1:-1, 1:-1, 2:] + fp[1:-1, 1:-1, :-2] - 6 * f) / h ** 2

    worst = 0.0
    for k in range(steps):
        on = ast > 0 and (k + 0.5) / steps >= ab_on
        idx = np.where(started)[0]
        Gn, res = (step1 if on else step0)(jnp.asarray(G[idx]), jnp.asarray(c.ravel()[idx]))
        G[idx] = np.asarray(Gn); worst = max(worst, float(jnp.max(res)))
        # phase
        Psi = np.zeros(phi.size)
        Psi[idx] = psi_pm(G[idx], A, ast if on else 0.0, b)
        Psi = Psi.reshape(phi.shape)
        force = beta * lap(phi) - 2 * kphi / (1 - phi) ** 3 - lam * 2 * phi * (Psi - psi_ref)
        phi = np.clip(phi + dt / eta_phi * force, 1e-6, 1 - 1e-6)
        newly = (~started) & (phi.ravel() > phi_min)
        G[newly] = g; started |= newly
        # nutrient and growth
        living = (G[:, 0:5] * G[:, 6:11]).sum(1).reshape(phi.shape) * started.reshape(phi.shape)
        rate = L0 * phi * living * np.maximum(Y * np.log(np.maximum(c, 1e-12)) - mu_b, 0.0)
        c = np.maximum(c + dt * (D * lap(c) - Y * rate), 1e-12); c[-1] = 1.0
        xi += dt * rate
    alpha = 1.0 + xi
    U = f3.elasticity(phi, alpha, n)
    return x, X, Z, phi, alpha, c, U, living, worst


def reference_energy(steps=80, T=0.2375, s=1.0):
    """Time average of the point model's stored energy over T* = 0..1 without antibiotic, from the
    calibration's start: the formation energy that keeps an untreated biofilm phase steady on
    average.  The absolute level of Psi_pm is not fixed by the calibration (only its derivatives
    enter the point model), so with lam > 0 a reference is needed; this is one choice."""
    g0, th = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1; g[6:11] = np.clip(g[6:11], 1e-10, 1 - 1e-10)
    A = np.asarray(ch.theta_to_matrices(jnp.asarray(th))[0])
    st = f3s.make_pm_step(th, s * T / steps, 0.0, np.ones(5))
    G = g[None].copy(); vals = []
    for _ in range(steps):
        G, _ = st(jnp.asarray(G), jnp.asarray([1.0])); G = np.asarray(G)
        vals.append(psi_pm(G, A, 0.0, np.ones(5))[0])
    return float(np.mean(vals))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    pref = reference_energy()
    print(f"reference energy (time average of Psi_pm, untreated) = {pref:.3f}")
    out = {}
    for lam, pr in ((0.0, 0.0), (1.0, 0.0), (1.0, pref)):
        for name, kw in (("no antibiotic", {}), ("antibiotic", dict(ast=0.65))):
            x, X, Z, phi, alpha, c, U, living, worst = run(n=n, lam=lam, psi_ref=pr, **kw)
            j = n // 2; h = 1.0 / n
            vol = phi.sum() * h ** 3
            d30 = [f3s.edge_disp(x, U, 30, az) for az in (0.0, np.pi)]
            grown = (phi * alpha ** 3).sum() * h ** 3
            print(f"lam {lam:.0f}, ref {pr:.2f}, {name:14s}: biofilm volume (reference) {vol:.5f}, grown volume {grown:.5f}, "
                  f"phi at centre {phi[j, j, 1]:.3f}, living at centre {living[j, j, 1]:.3f}, alpha-1 nutrient side "
                  f"{alpha[int(0.7 * n), j, 1] - 1:.3f}; edge 30 deg {d30[0]:.4f} / {d30[1]:.4f}; residual <= {worst:.1e}", flush=True)
            if not (lam > 0 and pr == 0.0):
                out[(lam, name)] = (phi[:, j, :], alpha[:, j, :])
    import figstyle
    figstyle.apply()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(8, 5.2), sharex=True, sharey=True)
    for a, (key, (ph, al)) in zip(ax.ravel(), out.items()):
        m = a.contourf(X[:, n // 2, :], Z[:, n // 2, :], ph, levels=np.linspace(0, 1, 11), cmap="Greens")
        a.contour(X[:, n // 2, :], Z[:, n // 2, :], al - 1, levels=[0.05, 0.1], colors="k", linewidths=0.8)
        a.set(title=(f"$\\lambda={key[0]:.0f}$, reference energy, " if key[0] > 0 else "$\\lambda=0$, ") + key[1], xlim=(0.1, 0.9), ylim=(0, 0.45), aspect="equal")
    for a in ax[1]:
        a.set_xlabel("$x$")
    for a in ax[:, 0]:
        a.set_ylabel("$z$")
    fig.colorbar(m, ax=ax, label=r"phase $\phi$ (lines: $\alpha-1=0.05,\ 0.1$)", shrink=0.8)
    fig.savefig(ROOT / "assets/fig_fuel_front_3d_coupled.pdf", bbox_inches="tight")
    print("wrote assets/fig_fuel_front_3d_coupled.pdf")
