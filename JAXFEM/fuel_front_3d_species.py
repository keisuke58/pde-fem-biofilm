"""fuel_front_3d_species.py -- fuel_front_3d.py with the calibrated five-species point model in place
of the lumped living share.

The composition is the two-level model in its decoupled limit (lam -> 0, beta_s = 0): at every node
inside the biofilm the point model of Klempt et al. 2026 runs on its own clock (s = 0.2375 per T*,
the calibration's time link), from the CH Day-1 state put to the inner liquid s_0 = 0.1.  The
living biomass is the catalyst of the growth reaction:
    xidot = L0 * phi * sum_i s_i psi_i * (Y ln c - mu_b),   cdot = D lap c - Y xidot,
mechanics as fuel_front_3d.py.  Each point-model step is solved with damped Newton to |Q| < 1e-10
(the residual of core_hamilton_1d, gate off; 60 damped Newton iterations per step).  Cases: no antibiotic; antibiotic from T* = 0.5 with
the same susceptibility for all species; antibiotic from T* = 0.5 acting on the anaerobes
(Vei, Fn, Pg ten times more susceptible).  The antibiotic parameters are assumed.
Output: assets/fig_fuel_front_3d_species.pdf.
"""
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)
ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "JAXFEM", ROOT / "ansys_usermat/coupling", ROOT / "ansys_usermat"):
    sys.path.insert(0, str(p))
import core_hamilton_1d as ch  # noqa: E402
import ecology_jax as ej  # noqa: E402
import fuel_front_3d as f3  # noqa: E402
import unified_1d as u1  # noqa: E402


def make_pm_step(theta, dt_pm, ast, b):
    A, _ = ch.theta_to_matrices(jnp.asarray(theta))
    p = jax.tree_util.tree_map(jnp.asarray, ej.default_hparams(dt_pm))
    p["A"] = A; p["b_diag"] = jnp.asarray(b, float); p["alpha"] = jnp.float64(ast)
    mask = jnp.ones(5, dtype=jnp.int64)

    def one(gp):
        F = lambda g: ch.residual(g, gp, p)

        def body(g, _):
            d = jnp.linalg.solve(jax.jacfwd(F)(g), -F(g))
            v, dv = g[0:11], d[0:11]
            lim = jnp.where(dv < 0, -0.9 * v / jnp.minimum(dv, -1e-300), jnp.where(dv > 0, 0.9 * (1 - v) / jnp.maximum(dv, 1e-300), jnp.inf))
            return g + jnp.minimum(1.0, jnp.min(lim)) * d, None

        g, _ = jax.lax.scan(body, ch.clip_state(gp, mask), jnp.arange(60))
        return g, jnp.max(jnp.abs(F(g)))
    return jax.jit(jax.vmap(one))


def run(n=24, D=0.2, L0=2.0, Y=1.0, c_stop=0.1, T=0.2375, s=0.2375, ast=0.0, ab_on=0.5, b=None):
    g0, th = u1.ch_case()
    g = g0.copy(); g[0:5] *= 0.9 / g[0:5].sum(); g[5] = 0.1; g[6:11] = np.clip(g[6:11], 1e-10, 1 - 1e-10)
    h = 1.0 / n
    x = np.linspace(0, 1, n + 1)
    X, Yg, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt((X - 0.5) ** 2 + (Yg - 0.5) ** 2 + Z ** 2)
    phi = 0.5 * (1 - np.tanh((r - 0.25) / 0.03)) * 0.9 + 1e-4
    act = phi > 0.01
    G = np.tile(g, (act.sum(), 1))
    dt = 0.15 * h * h / D; steps = int(np.ceil(T / dt)); dt = T / steps
    b = np.ones(5) if b is None else np.asarray(b, float)
    step0, step1 = make_pm_step(th, s * dt, 0.0, b), make_pm_step(th, s * dt, ast, b)
    c = np.full_like(phi, c_stop); c[-1] = 1.0; xi = np.zeros_like(phi); mu_b = Y * np.log(c_stop)
    worst = 0.0
    for k in range(steps):
        on = ast > 0 and (k + 0.5) / steps >= ab_on
        Gj, res = (step1 if on else step0)(jnp.asarray(G))
        G = np.asarray(Gj); worst = max(worst, float(jnp.max(res)))
        if k < 3:
            print(f'  step {k}: point-model residual {float(jnp.max(res)):.1e}', flush=True)
        living = np.zeros_like(phi); living[act] = (G[:, 0:5] * G[:, 6:11]).sum(1)
        rate = L0 * phi * living * np.maximum(Y * np.log(np.maximum(c, 1e-12)) - mu_b, 0.0)
        cp = np.pad(c, 1, mode="edge")
        lap = (cp[2:, 1:-1, 1:-1] + cp[:-2, 1:-1, 1:-1] + cp[1:-1, 2:, 1:-1] + cp[1:-1, :-2, 1:-1]
               + cp[1:-1, 1:-1, 2:] + cp[1:-1, 1:-1, :-2] - 6 * c) / h ** 2
        c = np.maximum(c + dt * (D * lap - Y * rate), 1e-12); c[-1] = 1.0
        xi += dt * rate
    alpha = 1.0 + xi
    U = f3.elasticity(phi, alpha, n)
    return x, X, Z, phi, alpha, c, U, living, worst


def edge_disp(x, U, el, az):
    from scipy.interpolate import RegularGridInterpolator
    itp = [RegularGridInterpolator((x, x, x), U[..., d]) for d in range(3)]
    e = np.radians(el)
    dvec = np.array([np.cos(e) * np.cos(az), np.cos(e) * np.sin(az), np.sin(e)])
    p = np.array([0.5, 0.5, 0.0]) + 0.25 * dvec
    return np.array([f(p)[0] for f in itp]) @ dvec


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    cases = (("no antibiotic", dict()),
             ("antibiotic, all species", dict(ast=0.65)),
             ("antibiotic, anaerobes", dict(ast=0.32, b=[1, 1, 10, 10, 10])))
    out = {}
    for name, kw in cases:
        x, X, Z, phi, alpha, c, U, living, worst = run(n=n, **kw)
        j = n // 2
        d30 = [edge_disp(x, U, 30, az) for az in (0.0, np.pi)]
        print(f"{name:24s}: living share at the colony centre {living[j, j, 1]:.3f}, alpha - 1 nutrient side "
              f"{alpha[int(0.7 * n), j, 1] - 1:.3f}, far side {alpha[int(0.3 * n), j, 1] - 1:.3f}; edge 30 deg: "
              f"nutrient side {d30[0]:.4f}, far side {d30[1]:.4f}; point-model residual <= {worst:.1e}", flush=True)
        out[name] = (alpha[:, j, :], U[:, j, :, :])
    import figstyle
    figstyle.apply()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.interpolate import RegularGridInterpolator
    th = np.linspace(0, np.pi, 91); P = np.stack([0.5 + 0.25 * np.cos(th), 0.25 * np.sin(th)], 1)
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0), sharey=True)
    vmax = max(o[0].max() for o in out.values()) - 1
    for a, (name, (al, Us)) in zip(ax, out.items()):
        a.contourf(X[:, n // 2, :], Z[:, n // 2, :], al - 1, 20, cmap="Oranges", vmin=0, vmax=vmax)
        iu = RegularGridInterpolator((x, x), Us[..., 0]); iw = RegularGridInterpolator((x, x), Us[..., 2])
        Q = P + np.stack([iu(P), iw(P)], 1)
        a.plot(P[:, 0], P[:, 1], "k--", lw=1); a.plot(Q[:, 0], Q[:, 1], "k-", lw=1.5)
        a.set(title=name, xlabel="$x$", aspect="equal", xlim=(0.15, 0.85), ylim=(0, 0.4))
    ax[0].set_ylabel("$z$ (section $y=0.5$)")
    fig.tight_layout(); fig.savefig(ROOT / "assets/fig_fuel_front_3d_species.pdf")
    print("wrote assets/fig_fuel_front_3d_species.pdf")
