"""
benchmark_klempt2024.py
=======================
Phase 0: Klempt 2024 単一菌種ベンチマーク

Klempt et al. (2024) Biomech Model Mechanobiol, Case 1:
  - Domain  : 20×20 μm (normalized to [0,1]×[0,1])
  - Biofilm : 中心に小さなシード (φ₁ = 0.3 in 5×5 patch)
  - Nutrient: 右上コーナー Dirichlet c=1, 他 Neumann
  - Species : S.oralis のみ active (active_mask=[1,0,0,0,0])
  - η = 0 (粘性なし → F = Fe·Fg, Klemptと同じ)

合格基準: φ場の最終形態がKlempt Fig 2と定性一致
  (バイオフィルムが栄養素源コーナーに向かって指向性成長)

Run:
    python benchmark_klempt2024.py
    python benchmark_klempt2024.py --nx 40 --n_macro 200 --save
"""

import argparse
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from core_hamilton_2d_nutrient import (
    Config2D,
    _make_reaction_step,
    _make_nutrient_step_mixed,
    diffusion_step_species_2d,
    clip_state,
)

# ---------------------------------------------------------------------------
# Klempt 2024 parameters (Table 1 / normalized)
# ---------------------------------------------------------------------------
KLEMPT_PARAMS = {
    "Kp1"        : 1e-4,   # penalty
    "c_hamilton" : 100.0,  # growth driving force scale
    "alpha"      : 100.0,  # Hamilton alpha constant
    "K_hill"     : 0.0,    # Hill gate off (single species)
    "n_hill"     : 4.0,
    # Diffusion — low D_c creates strong gradient (Klempt: diffusion limited)
    "D_eff_s1"   : 2e-4,   # S.oralis (slow: stays localized)
    "D_c"        : 5e-4,   # nutrient diffusion (small → gradient forms)
    "k_monod"    : 0.5,    # Monod half-saturation
    "g_s1"       : 8.0,    # consumption rate (high → depletion near biofilm)
    # Single-species interaction
    "a11"        : 2.0,    # self-competition
    "b1"         : 1.5,    # growth rate
}


def make_initial_state_klempt(Nx, Ny, seed_radius=0.08):
    """
    初期条件: Klempt Case 1
      φ₁ = 0.3 in central seed patch, 0 elsewhere
      φ₀ = 1 - φ₁ (void)
      c  = 1.0 uniform (栄養素は初期均一)
    """
    x = jnp.linspace(0, 1, Nx)
    y = jnp.linspace(0, 1, Ny)
    xv, yv = jnp.meshgrid(x, y, indexing="ij")

    # Circular seed at center
    r = jnp.sqrt((xv - 0.5)**2 + (yv - 0.5)**2)
    phi1 = jnp.where(r < seed_radius, 0.3, 0.0)

    G_2d = jnp.zeros((Nx, Ny, 12), dtype=jnp.float64)
    G_2d = G_2d.at[:, :, 0].set(phi1)           # S.oralis
    G_2d = G_2d.at[:, :, 5].set(1.0 - phi1)     # void φ₀
    G_2d = G_2d.at[:, :, 6].set(                 # ψ₁
        jnp.where(phi1 > 1e-6, 0.999, 0.0)
    )
    G_2d = G_2d.at[:, :, 11].set(0.0)           # γ

    G_flat = G_2d.reshape(Nx * Ny, 12)
    c = jnp.ones((Nx, Ny), dtype=jnp.float64)
    return G_flat, c


def build_single_species_params(cfg, p=KLEMPT_PARAMS):
    """active_mask=[1,0,0,0,0], A=[[a11,0,...]], b=[b1,0,...]"""
    active_mask = jnp.array([1, 0, 0, 0, 0], dtype=jnp.int64)
    A = jnp.zeros((5, 5), dtype=jnp.float64).at[0, 0].set(p["a11"])
    b_diag = jnp.zeros(5, dtype=jnp.float64).at[0].set(p["b1"])
    return {
        "dt_h"       : cfg.dt_h,
        "Kp1"        : cfg.Kp1,
        "Eta"        : jnp.ones(5, dtype=jnp.float64),
        "EtaPhi"     : jnp.ones(5, dtype=jnp.float64),
        "c"          : cfg.c_hamilton,
        "alpha"      : cfg.alpha,
        "K_hill"     : jnp.array(0.0),
        "n_hill"     : jnp.array(4.0),
        "A"          : A,
        "b_diag"     : b_diag,
        "active_mask": active_mask,
    }


def run_klempt_benchmark(Nx=30, Ny=30, n_macro=150, save=False):
    p = KLEMPT_PARAMS

    cfg = Config2D(
        Nx=Nx, Ny=Ny,
        Lx=1.0, Ly=1.0,
        dt_h=1e-5,
        n_react_sub=20,
        n_macro=n_macro,
        save_every=max(1, n_macro // 6),
        D_eff=np.array([p["D_eff_s1"], 0, 0, 0, 0]),
        D_c=p["D_c"],
        k_monod=p["k_monod"],
        g_consumption=np.array([p["g_s1"], 0, 0, 0, 0]),
        c_boundary=1.0,
        Kp1=p["Kp1"],
        c_hamilton=p["c_hamilton"],
        alpha=p["alpha"],
        K_hill=p["K_hill"],
        n_hill=p["n_hill"],
        newton_iters=6,
    )

    params = build_single_species_params(cfg, p)
    _reaction_step = _make_reaction_step(cfg.n_react_sub, cfg.newton_iters)
    _nutrient_step = _make_nutrient_step_mixed(n_sub_c=5)

    G, c = make_initial_state_klempt(Nx, Ny)

    phi_snaps, c_snaps, t_snaps = [np.array(G.reshape(Nx, Ny, 12)[:, :, 0])], [np.array(c)], [0.0]

    print(f"\n{'='*55}")
    print(f"Klempt 2024 Phase 0 Benchmark | Nx={Nx} Ny={Ny}")
    print(f"  single species (S.oralis), eta=0, a11={p['a11']}, b1={p['b1']}")
    print(f"  n_macro={n_macro}, dt_macro={cfg.dt_macro:.2e}")
    print(f"{'='*55}\n")

    for step in range(1, n_macro + 1):
        # (1) Hamilton reaction
        G = _reaction_step(G, params)

        # (2) Species diffusion (only species 0 active)
        # diffusion_step_species_2d expects (Nx, Ny, 5)
        phi_2d = G.reshape(Nx, Ny, 12)[:, :, :5]   # (Nx, Ny, 5)
        D_eff = jnp.array(cfg.D_eff)
        phi_2d = diffusion_step_species_2d(phi_2d, D_eff, cfg.dt_macro, cfg.dx, cfg.dy)

        # Write back diffused phi
        G = G.reshape(Nx, Ny, 12).at[:, :, :5].set(phi_2d).reshape(Nx * Ny, 12)

        # _make_nutrient_step_mixed returns step(c, phi_2d, D_c, k_M, g_cons, c_bc, dx, dy, dt)
        # phi_2d is (Nx, Ny, 5)
        c = _nutrient_step(
            c, phi_2d,
            cfg.D_c, cfg.k_monod,
            np.array(cfg.g_consumption), cfg.c_boundary,
            cfg.dx, cfg.dy, cfg.dt_macro,
        )

        # Enforce corner Dirichlet: top-right = 1.0
        c = c.at[-1, -1].set(1.0)

        if step % cfg.save_every == 0:
            phi1 = np.array(G.reshape(Nx, Ny, 12)[:, :, 0])
            print(f"  step {step:4d}/{n_macro}  φ_max={phi1.max():.3f}  c_min={float(c.min()):.3f}")
            phi_snaps.append(phi1)
            c_snaps.append(np.array(c))
            t_snaps.append(step * cfg.dt_macro)

    _plot_results(phi_snaps, c_snaps, t_snaps, Nx, Ny, save)
    return phi_snaps, c_snaps, t_snaps


def _plot_results(phi_snaps, c_snaps, t_snaps, Nx, Ny, save):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    n = len(phi_snaps)
    fig = plt.figure(figsize=(4 * n, 8))
    gs = gridspec.GridSpec(2, n, hspace=0.4, wspace=0.3)

    for i, (phi, c, t) in enumerate(zip(phi_snaps, c_snaps, t_snaps)):
        ax1 = fig.add_subplot(gs[0, i])
        im1 = ax1.imshow(phi.T, origin="lower", cmap="YlOrBr", vmin=0, vmax=0.5)
        ax1.set_title(f"φ₁  t={t:.3f}")
        ax1.set_xlabel("x"); ax1.set_ylabel("y")
        plt.colorbar(im1, ax=ax1, fraction=0.046)

        ax2 = fig.add_subplot(gs[1, i])
        im2 = ax2.imshow(c.T, origin="lower", cmap="Blues", vmin=0, vmax=1)
        ax2.set_title(f"c  t={t:.3f}")
        ax2.set_xlabel("x"); ax2.set_ylabel("y")
        plt.colorbar(im2, ax=ax2, fraction=0.046)

    fig.suptitle("Klempt 2024 Phase 0 Benchmark: Single-species directional growth\n"
                 "(Nutrient source: top-right corner → biofilm grows toward it)",
                 fontsize=11)

    if save:
        out = Path(__file__).parent / "benchmark_klempt2024_results.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"\nSaved: {out}")
    else:
        plt.show()

    plt.close(fig)

    # Pass/Fail check
    # BC: Dirichlet c=1 at top edge (y=Ly) → nutrient comes from top
    # Expect: φ mean in top half > φ mean in bottom half
    phi_final = phi_snaps[-1]
    Nx, Ny = phi_final.shape
    phi_top    = phi_final[:, Ny//2:].mean()   # high-nutrient side
    phi_bottom = phi_final[:, :Ny//2].mean()   # low-nutrient side
    c_final = c_snaps[-1]
    c_top    = c_final[:, Ny//2:].mean()
    c_bottom = c_final[:, :Ny//2].mean()
    passed = phi_top > phi_bottom and c_top > c_bottom
    print(f"\n{'='*55}")
    print(f"BENCHMARK RESULT: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"  φ_mean top (nutrient side)  : {phi_top:.4f}")
    print(f"  φ_mean bottom               : {phi_bottom:.4f}")
    print(f"  c_mean top                  : {c_top:.4f}")
    print(f"  c_mean bottom               : {c_bottom:.4f}")
    print(f"  → Biofilm enriched near nutrient source: {passed}")
    print(f"\n  NOTE: NSP dynamics ≠ Klempt α-growth PDE.")
    print(f"  True Klempt match requires implementing dα/dt equation.")
    print(f"  Phase 0 verifies: nutrient modulation of NSP is directional.")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Klempt 2024 Phase 0 Benchmark")
    parser.add_argument("--nx", type=int, default=30)
    parser.add_argument("--ny", type=int, default=30)
    parser.add_argument("--n_macro", type=int, default=150)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()

    run_klempt_benchmark(Nx=args.nx, Ny=args.ny, n_macro=args.n_macro, save=args.save)
