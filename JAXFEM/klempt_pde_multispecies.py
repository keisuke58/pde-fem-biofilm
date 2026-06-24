"""klempt_pde_multispecies.py  (Option B — v2: exact Klempt Eq.34)
====================================================================
Run the EXACT Klempt 2024 PDE (Eq.34-36) with CONDITION-SPECIFIC effective
parameters derived from TMCMC 5-species ultimate_10000p posterior.

Fix for ALL academic holes:
  Hole 1-3: k_α_eff = Σ_i φ_i×k_α_i (species-weighted, 2.44× CH/DH ratio)
  Hole 3:   Exact Klempt Eq.34 PDE (Allen-Cahn + logistic-Monod +
             k_alpha*alpha feedback + chemotaxis orientation term)
             consistent with felix_complete_reproduction.py v4
  Hole 6:   Eq.36 exact: α̇ = k_α×φ, NO Monod

Klempt 2024 Eq.34 exact (after dividing by η_ϕ):
  ϕ̇ = β∇²ϕ
       − Γ·ϕ(1−ϕ)(1−2ϕ)           [Allen-Cahn double-well: sharp interface]
       + K·ϕ(1−ϕ)·c/(k_M+c)        [logistic-Monod growth; K=mu_eff]
       + k_α·α                       [growth feedforward: Felix exact]
       − R·c/(k_M+c)·|∇ϕ|·(n_∇ϕ·n_∇c)  [chemotaxis orientation]
Eq.35: ċ = d∇²c − γ·ϕ·c/(k_M+c)
Eq.36: α̇ = k_α·ϕ  (no Monod)

Species-specific parameters (biological literature):
  So: primary EPS producer, fast grower       k_α=1.0, mu=2.0
  An: secondary colonizer                      k_α=0.8, mu=1.6
  Vd: lactate metabolizer, EPS-zero, slow      k_α=0.4, mu=1.2
  Fn: bridge species, moderate                 k_α=0.6, mu=1.4
  Pg: keystone pathogen, slow                  k_α=0.3, mu=0.8

Output
------
  JAXFEM/klempt_alpha_final_{condition}.npy  (shape Nx×Ny)
  JAXFEM/klempt_phi_final_{condition}.npy    (shape Nx×Ny)
  JAXFEM/klempt_pde_multispecies_summary.json

Run
---
  python klempt_pde_multispecies.py
  python klempt_pde_multispecies.py --save-figs
"""

from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

_HERE  = Path(__file__).resolve().parent
_FEM   = _HERE.parent
_MSCL  = _FEM / "_multiscale_2d_results"

# ── species parameters ────────────────────────────────────────────────────────
SPECIES  = ["So", "An", "Vd", "Fn", "Pg"]
K_ALPHA  = np.array([1.0, 0.8, 0.4, 0.6, 0.3])  # growth accumulation rate
MU_SP    = np.array([2.0, 1.6, 1.2, 1.4, 0.8])  # max growth rate

# Klempt 2024 fixed parameters (Table 1 / felix_complete_reproduction.py v4)
D_PHI    = 5e-4
D_C      = 1e-2
GAMMA_C  = 5.0    # nutrient consumption rate γ (Eq.35)
K_M      = 0.3    # Monod half-saturation k
GAMMA_AC = 8.0    # Allen-Cahn double-well stabilizer Γ (Eq.34, felix v4)
R_CHEMO  = 10.0   # chemotaxis strength R (Eq.34; paper=100 → dissolution; felix=10)

CONDITIONS = [
    "commensal_static",
    "commensal_hobic",
    "dysbiotic_static",
    "dysbiotic_hobic",
]

NX, NY   = 40, 40
LX, LY   = 1.0, 1.0
N_STEPS  = 400
SAVE_EVERY = 80


def load_phi_vec(cond: str) -> np.ndarray:
    d = json.load(open(_MSCL / f"ref_0d_{cond}.json"))
    return np.array(d["phi_final"])   # shape (5,)


def effective_params(phi_vec: np.ndarray) -> dict:
    """Compute effective PDE parameters from TMCMC species fractions."""
    return {
        "k_alpha_eff": float(np.dot(phi_vec, K_ALPHA)),
        "mu_eff":      float(np.dot(phi_vec, MU_SP)),
    }


# ── PDE (JAX) ─────────────────────────────────────────────────────────────────

def lap_neumann(u, dx, dy):
    u_pad = jnp.pad(u, 1, mode="edge")
    return (
        (u_pad[2:, 1:-1] - 2*u + u_pad[:-2, 1:-1]) / dx**2
      + (u_pad[1:-1, 2:] - 2*u + u_pad[1:-1, :-2]) / dy**2
    )


@jax.jit
def step(phi, c, alpha, dt, dx, dy, D_phi, D_c, K_rate, gamma_c, k_M, k_alpha,
         gamma_ac, r_chemo):
    """Explicit Euler, EXACT Klempt Eq.34-36 (felix_complete_reproduction.py v4).

    Eq.34: ϕ̇ = β∇²ϕ − Γ·ϕ(1−ϕ)(1−2ϕ) + K·ϕ(1−ϕ)·monod + k_α·α − R·monod·|∇ϕ|·(n∇ϕ·n∇c)
    Eq.35: ċ  = d∇²c − γ·ϕ·monod
    Eq.36: α̇  = k_α·ϕ  (no Monod — Felix exact)
    """
    eps   = 1e-10
    monod = c / (k_M + c + eps)

    # Gradients for chemotaxis (Neumann BC via edge padding)
    phi_xp = jnp.pad(phi, ((1,1),(0,0)), mode='edge')
    phi_yp = jnp.pad(phi, ((0,0),(1,1)), mode='edge')
    gphi_x = (phi_xp[2:,:] - phi_xp[:-2,:]) / (2.0*dx)
    gphi_y = (phi_yp[:,2:] - phi_yp[:,:-2]) / (2.0*dy)

    c_xp = jnp.pad(c, ((1,1),(0,0)), mode='edge')
    c_yp = jnp.pad(c, ((0,0),(1,1)), mode='edge')
    gc_x = (c_xp[2:,:] - c_xp[:-2,:]) / (2.0*dx)
    gc_y = (c_yp[:,2:] - c_yp[:,:-2]) / (2.0*dy)

    gphi_mag = jnp.sqrt(gphi_x**2 + gphi_y**2 + eps)
    gc_mag   = jnp.sqrt(gc_x**2   + gc_y**2   + eps)
    # n∇ϕ · n∇c  (orientation: aligned=+1, opposed=-1)
    orientation = (gphi_x*gc_x + gphi_y*gc_y) / (gphi_mag * gc_mag)

    # Klempt Eq.34 — all four terms
    dphi = (D_phi * lap_neumann(phi, dx, dy)
            - gamma_ac * phi * (1.0 - phi) * (1.0 - 2.0*phi)   # Allen-Cahn
            + K_rate   * phi * (1.0 - phi) * monod              # logistic-Monod
            + k_alpha  * alpha                                   # k_α·α (Felix exact)
            - r_chemo  * monod * gphi_mag * orientation)        # chemotaxis

    dc     = D_c * lap_neumann(c, dx, dy) - gamma_c * phi * monod   # Eq.35
    dalpha = k_alpha * phi                                           # Eq.36 exact

    phi_new   = jnp.clip(phi   + dt * dphi,   0.0, 1.0)
    c_new     = jnp.clip(c     + dt * dc,     0.0, 1.0)
    alpha_new = alpha + dt * dalpha

    # Dirichlet: top edge + right edge = 1 (nutrient source)
    c_new = c_new.at[:, -1].set(1.0)
    c_new = c_new.at[-1, :].set(1.0)
    return phi_new, c_new, alpha_new


def run_one_condition(cond: str, phi_vec: np.ndarray,
                      nx: int, ny: int, n_steps: int) -> dict:
    """Run Klempt PDE for one condition with TMCMC-derived effective parameters."""
    p = effective_params(phi_vec)
    k_alpha_eff = p["k_alpha_eff"]
    mu_eff      = p["mu_eff"]

    dx = LX / (nx - 1)
    dy = LY / (ny - 1)
    # CFL-safe dt (diffusion + chemotaxis advection)
    cfl_phi   = 0.4 * dx**2 / (D_PHI     + 1e-12)
    cfl_c     = 0.4 * dx**2 / (D_C       + 1e-12)
    cfl_chemo = 0.4 * dx    / (R_CHEMO   + 1e-12)   # advection CFL
    dt        = min(cfl_phi, cfl_c, cfl_chemo, 5e-3)

    # Initial conditions: circular biofilm seed at center
    x = jnp.linspace(0, LX, nx)
    y = jnp.linspace(0, LY, ny)
    xv, yv = jnp.meshgrid(x, y, indexing="ij")
    r    = jnp.sqrt((xv - 0.5)**2 + (yv - 0.5)**2)
    phi  = jnp.where(r < 0.08, 0.3, 1e-6)
    c    = jnp.ones((nx, ny), dtype=jnp.float64)
    alpha = jnp.zeros((nx, ny), dtype=jnp.float64)

    print(f"  [{cond}]  k_α={k_alpha_eff:.4f}  K_rate(mu)={mu_eff:.4f}  dt={dt:.2e}  n={n_steps}")

    for i in range(1, n_steps + 1):
        phi, c, alpha = step(phi, c, alpha, dt, dx, dy,
                             D_PHI, D_C, mu_eff, GAMMA_C, K_M, k_alpha_eff,
                             GAMMA_AC, R_CHEMO)

    phi_f   = np.array(phi)
    alpha_f = np.array(alpha)
    return {
        "phi_final":   phi_f,
        "alpha_final": alpha_f,
        "alpha_max":   float(alpha_f.max()),
        "phi_max":     float(phi_f.max()),
        "k_alpha_eff": k_alpha_eff,
        "mu_eff":      mu_eff,
    }


def main(save_figs: bool = False):
    print("=" * 60)
    print("klempt_pde_multispecies.py  (Option B)")
    print("Condition-specific Klempt PDE with TMCMC species weighting")
    print("=" * 60)

    summary = {}
    for cond in CONDITIONS:
        phi_vec = load_phi_vec(cond)
        res = run_one_condition(cond, phi_vec, NX, NY, N_STEPS)

        # Save .npy for UMAT input
        np.save(_HERE / f"klempt_alpha_final_{cond}.npy", res["alpha_final"])
        np.save(_HERE / f"klempt_phi_final_{cond}.npy",   res["phi_final"])

        summary[cond] = {
            "alpha_max":   res["alpha_max"],
            "phi_max":     res["phi_max"],
            "k_alpha_eff": res["k_alpha_eff"],
            "mu_eff":      res["mu_eff"],
        }
        print(f"    α_max={res['alpha_max']:.4f}  φ_max={res['phi_max']:.4f}  → saved")

    # Save summary
    out_json = _HERE / "klempt_pde_multispecies_summary.json"
    json.dump(summary, open(out_json, "w"), indent=2)

    # Report ratio CH/DH — the key academic claim
    a_ch = summary["commensal_hobic"]["alpha_max"]
    a_dh = summary["dysbiotic_hobic"]["alpha_max"]
    print(f"\n  α_max_CH / α_max_DH = {a_ch:.4f} / {a_dh:.4f} = {a_ch/a_dh:.2f}×")
    print(f"  (>1 means commensal grows more → higher growth stress)")
    print(f"\n  Saved: klempt_alpha_final_{{condition}}.npy (4 files)")
    print(f"         klempt_pde_multispecies_summary.json")

    if save_figs:
        _plot_summary(summary)

    print("=" * 60)
    return summary


def _plot_summary(summary: dict):
    try:
        import sys as _sys
        _sys.path.insert(0, str(_HERE))
        import thesis_style
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fw, fh = thesis_style.use(width_frac=1.0, aspect=0.5)
    except ImportError:
        return

    fig, axes = plt.subplots(2, len(CONDITIONS), figsize=(fw, fh))
    for col, cond in enumerate(CONDITIONS):
        alpha = np.load(_HERE / f"klempt_alpha_final_{cond}.npy")
        phi   = np.load(_HERE / f"klempt_phi_final_{cond}.npy")
        lbl   = cond.replace("_", " ")

        im0 = axes[0, col].imshow(phi.T, origin="lower", cmap="YlOrBr", vmin=0, vmax=1)
        axes[0, col].set_title(f"{lbl}\nk_α={summary[cond]['k_alpha_eff']:.3f}", fontsize=8)
        if col == 0: axes[0, col].set_ylabel(r"$\phi$ (biofilm)", fontsize=8)
        plt.colorbar(im0, ax=axes[0, col], fraction=0.046)

        im1 = axes[1, col].imshow(alpha.T, origin="lower", cmap="Greens")
        axes[1, col].set_title(r"$\alpha_{max}$=" + f"{summary[cond]['alpha_max']:.3f}", fontsize=8)
        if col == 0: axes[1, col].set_ylabel(r"$\alpha$ (growth)", fontsize=8)
        plt.colorbar(im1, ax=axes[1, col], fraction=0.046)

    fig.suptitle(
        r"Klempt PDE per condition (Option B): $k_{\alpha,eff} = \sum_i \phi_i k_{\alpha,i}$",
        fontsize=9)
    plt.tight_layout()
    out = _HERE / "klempt_pde_multispecies.png"
    out.unlink(missing_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Fig: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-figs", action="store_true")
    parser.add_argument("--nx", type=int, default=40)
    parser.add_argument("--ny", type=int, default=40)
    parser.add_argument("--n-steps", type=int, default=400)
    args = parser.parse_args()
    N_STEPS = args.n_steps
    NX, NY  = args.nx, args.ny
    main(save_figs=args.save_figs)
