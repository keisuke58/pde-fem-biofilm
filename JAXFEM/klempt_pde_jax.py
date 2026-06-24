"""
klempt_pde_jax.py
=================
Phase 0a: Klempt 2024 の支配方程式を JAX で直接実装

Klempt et al. (2024) Biomech Model Mechanobiol
https://pmc.ncbi.nlm.nih.gov/articles/PMC11554842/

State variables (single species):
  phi(x,y,t) : biofilm density  [0, 1]
  c(x,y,t)   : nutrient concentration [0, 1]
  alpha(x,y,t): local expansion parameter (cumulative growth)

Governing equations (simplified from Hamilton principle):
  dphi/dt = D_phi * lap(phi) + mu * phi * (1 - phi) * c / (k_M + c)
  dc/dt   = D_c  * lap(c)   - gamma * phi * c / (k_M + c)
  dalpha/dt = k_alpha * phi * c / (k_M + c)

Boundary conditions (Klempt Case 1):
  phi : Neumann (zero-flux) all walls
  c   : Dirichlet c=1 at top-right corner; Neumann elsewhere
  alpha: Neumann all walls

Run:
    python klempt_pde_jax.py
    python klempt_pde_jax.py --nx 40 --ny 40 --n_steps 500 --save
"""

import argparse
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


# ---------------------------------------------------------------------------
# Parameters (from Klempt 2024, normalized units)
# ---------------------------------------------------------------------------
class KlemptConfig:
    def __init__(self, Nx=40, Ny=40, Lx=1.0, Ly=1.0, n_steps=400, save_every=80):
        self.Nx = Nx
        self.Ny = Ny
        self.Lx = Lx
        self.Ly = Ly
        self.dx = Lx / (Nx - 1)
        self.dy = Ly / (Ny - 1)
        self.n_steps = n_steps
        self.save_every = save_every

        # PDE parameters
        self.D_phi   = 5e-4   # biofilm diffusion (slow: stays localized)
        self.D_c     = 1e-2   # nutrient diffusion
        self.mu      = 2.0    # max growth rate
        self.gamma   = 5.0    # nutrient consumption rate
        self.k_M     = 0.3    # Monod half-saturation
        self.k_alpha = 1.0    # expansion accumulation rate

        # CFL-safe dt
        cfl_phi = 0.4 * self.dx**2 / (self.D_phi + 1e-12)
        cfl_c   = 0.4 * self.dx**2 / (self.D_c   + 1e-12)
        self.dt = min(cfl_phi, cfl_c, 5e-3)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def lap_neumann(u, dx, dy):
    """2D Laplacian with Neumann (zero-flux) BCs via padding."""
    u_pad = jnp.pad(u, 1, mode="edge")
    return (
        (u_pad[2:, 1:-1] - 2 * u[...] + u_pad[:-2, 1:-1]) / dx**2
        + (u_pad[1:-1, 2:] - 2 * u[...] + u_pad[1:-1, :-2]) / dy**2
    )


def apply_corner_dirichlet(c, val=1.0):
    """Dirichlet c=val at top-right corner (Klempt Case 1 nutrient source).
    Apply entire top edge + right edge to create strong gradient."""
    c = c.at[:, -1].set(val)   # top edge (y=Ly): full Dirichlet
    c = c.at[-1, :].set(val)   # right edge (x=Lx): full Dirichlet
    return c


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

def make_initial(cfg, seed_r=0.08):
    """
    phi: circular seed at center (r < seed_r)
    c  : uniform = 1 everywhere initially
    alpha: 0 everywhere
    """
    x = jnp.linspace(0, cfg.Lx, cfg.Nx)
    y = jnp.linspace(0, cfg.Ly, cfg.Ny)
    xv, yv = jnp.meshgrid(x, y, indexing="ij")
    r = jnp.sqrt((xv - 0.5)**2 + (yv - 0.5)**2)
    phi   = jnp.where(r < seed_r, 0.3, 1e-6)
    c     = jnp.ones((cfg.Nx, cfg.Ny), dtype=jnp.float64)
    alpha = jnp.zeros((cfg.Nx, cfg.Ny), dtype=jnp.float64)
    return phi, c, alpha


# ---------------------------------------------------------------------------
# Single time step
# ---------------------------------------------------------------------------

@jax.jit
def step(phi, c, alpha, cfg_dt, cfg_dx, cfg_dy,
         D_phi, D_c, mu, gamma, k_M, k_alpha):
    """Explicit Euler step for Klempt φ-c-α PDE."""
    eps = 1e-10
    monod = c / (k_M + c + eps)

    # Laplacians (Neumann BCs)
    lap_phi = lap_neumann(phi, cfg_dx, cfg_dy)
    lap_c   = lap_neumann(c,   cfg_dx, cfg_dy)

    # RHS
    dphi   = D_phi * lap_phi + mu    * phi * (1.0 - phi) * monod
    dc     = D_c   * lap_c   - gamma * phi * monod
    dalpha = k_alpha * phi * monod

    phi_new   = jnp.clip(phi   + cfg_dt * dphi,   0.0, 1.0)
    c_new     = jnp.clip(c     + cfg_dt * dc,      0.0, 1.0)
    alpha_new = alpha + cfg_dt * dalpha

    # Enforce Dirichlet: top edge + right edge = 1 (Klempt Case 1 nutrient source)
    c_new = c_new.at[:, -1].set(1.0)
    c_new = c_new.at[-1, :].set(1.0)

    return phi_new, c_new, alpha_new


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------

def run(Nx=40, Ny=40, n_steps=400, save=False):
    cfg = KlemptConfig(Nx=Nx, Ny=Ny, n_steps=n_steps,
                       save_every=max(1, n_steps // 6))

    phi, c, alpha = make_initial(cfg)

    phi_snaps  = [np.array(phi)]
    c_snaps    = [np.array(c)]
    alpha_snaps = [np.array(alpha)]
    t_snaps    = [0.0]
    t = 0.0

    print(f"\n{'='*55}")
    print(f"Klempt PDE Phase 0a | Nx={Nx} Ny={Ny} dt={cfg.dt:.2e}")
    print(f"  D_phi={cfg.D_phi}  D_c={cfg.D_c}  mu={cfg.mu}  gamma={cfg.gamma}")
    print(f"{'='*55}\n")

    for i in range(1, n_steps + 1):
        phi, c, alpha = step(
            phi, c, alpha,
            cfg.dt, cfg.dx, cfg.dy,
            cfg.D_phi, cfg.D_c, cfg.mu, cfg.gamma, cfg.k_M, cfg.k_alpha,
        )
        t += cfg.dt

        if i % cfg.save_every == 0:
            print(f"  step {i:4d}/{n_steps}  φ_max={float(phi.max()):.3f}"
                  f"  c_min={float(c.min()):.3f}  α_max={float(alpha.max()):.3f}")
            phi_snaps.append(np.array(phi))
            c_snaps.append(np.array(c))
            alpha_snaps.append(np.array(alpha))
            t_snaps.append(float(t))

    _check_and_plot(phi_snaps, c_snaps, alpha_snaps, t_snaps, Nx, Ny, save)

    # Save final fields as .npy for Phase 1 stress analysis
    out_dir = Path(__file__).parent
    np.save(out_dir / "klempt_phi_final.npy",   phi_snaps[-1])
    np.save(out_dir / "klempt_c_final.npy",     c_snaps[-1])
    np.save(out_dir / "klempt_alpha_final.npy", alpha_snaps[-1])
    print(f"Saved npy: klempt_{{phi,c,alpha}}_final.npy  (shape {phi_snaps[-1].shape})")

    return phi_snaps, c_snaps, alpha_snaps


# ---------------------------------------------------------------------------
# Pass/Fail + plotting
# ---------------------------------------------------------------------------

def _check_and_plot(phi_snaps, c_snaps, alpha_snaps, t_snaps, Nx, Ny, save):
    phi_f = phi_snaps[-1]
    c_f   = c_snaps[-1]
    alpha_f = alpha_snaps[-1]

    # Nutrient source is at top-right corner (-1, -1)
    # Expect: phi and alpha accumulate in top-right quadrant
    q_tr = phi_f[Nx//2:, Ny//2:].mean()
    q_bl = phi_f[:Nx//2, :Ny//2].mean()
    a_tr = alpha_f[Nx//2:, Ny//2:].mean()
    a_bl = alpha_f[:Nx//2, :Ny//2].mean()
    c_min = c_f.min()

    passed = (q_tr > q_bl) and (a_tr > a_bl) and (c_min < 0.9)

    print(f"\n{'='*55}")
    print(f"BENCHMARK Phase 0a: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"  φ_mean top-right: {q_tr:.4f}  bottom-left: {q_bl:.4f}")
    print(f"  α_mean top-right: {a_tr:.4f}  bottom-left: {a_bl:.4f}")
    print(f"  c_min (depletion): {c_min:.4f}  [need < 0.9]")
    print(f"{'='*55}\n")

    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import thesis_style; fw, fh = thesis_style.use(width_frac=1.0, aspect=1.5)
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("matplotlib not available")
        return

    n = len(phi_snaps)
    fig = plt.figure(figsize=(fw, fh))
    gs = gridspec.GridSpec(3, n, hspace=0.45, wspace=0.3)

    for i, (phi, c, alpha, t) in enumerate(zip(phi_snaps, c_snaps, alpha_snaps, t_snaps)):
        for row, (data, cmap, label, vmax) in enumerate([
            (phi,   "YlOrBr", r"$\phi$ (biofilm)",    0.5),
            (c,     "Blues",  r"$c$ (nutrient)",       1.0),
            (alpha, "Greens", r"$\alpha$ (expansion)", None),
        ]):
            ax = fig.add_subplot(gs[row, i])
            im = ax.imshow(data.T, origin="lower", cmap=cmap,
                           vmin=0, vmax=vmax)
            ax.set_title(f"{label}\n$t={t:.3f}$")
            ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
            ax.plot([Nx//2], [Ny//2], "r+", ms=8)
            ax.plot([Nx-1], [Ny-1], "g*", ms=10)
            plt.colorbar(im, ax=ax, fraction=0.046)

    status = "PASS" if passed else "FAIL"
    fig.suptitle(
        r"Klempt 2024 Phase 0a: $\phi$-$c$-$\alpha$ PDE (single species) --- " + status + "\n"
        r"$\bigstar$ = nutrient corner (top-right),  $+$ = biofilm seed (center)",
        fontsize=9)

    if save:
        out = Path(__file__).parent / "klempt_pde_phase0a.png"
        out.unlink(missing_ok=True)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx",      type=int,  default=40)
    parser.add_argument("--ny",      type=int,  default=40)
    parser.add_argument("--n_steps", type=int,  default=400)
    parser.add_argument("--save",    action="store_true")
    args = parser.parse_args()
    run(Nx=args.nx, Ny=args.ny, n_steps=args.n_steps, save=args.save)
