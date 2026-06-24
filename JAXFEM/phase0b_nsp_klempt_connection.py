"""
phase0b_nsp_klempt_connection.py
=================================
Phase 0b: NSP Hamilton ODE ↔ Klempt 2024 PDE 接続確認

NSP ODE（菌種競合）と Klempt PDE（バイオフィルム密度）は
別々の方程式を解くが、マルチスケールモデルとして接続される：

  マクロスケール (Klempt PDE)
    ∂φ_total/∂t = D_φ∇²φ + μφ(1-φ)c/(k_M+c)
    ∂c/∂t = D_c∇²c - γφc/(k_M+c)
    dα/dt = k_α·φ·c/(k_M+c)
    → φ_total(x,y,t), c(x,y,t), α(x,y,t)

  メゾスケール (NSP ODE, single species → So のみ)
    active_mask = [1,0,0,0,0]   (S.oralis のみ)
    η = 0                         (粘性なし → F = Fe·Fg)
    c_local(t) = c_klempt(x,y,t)  (局所栄養素を Klempt から入力)
    → φ_So(t) ← 0D 競合ダイナミクス

  接続:
    φ_So_eq(x,y) = φ_total(x,y,T) × (So fraction from NSP)
    α(x,y) → ε_growth → UMAT/FEM 応力 (Phase 1 と同じ)

合格基準:
  1. NSP single-species (So only) が c_local → 正の平衡 φ_So > 0
  2. c を Klempt PDE から受け取ると φ_So が変化すること（栄養素依存性）

Run:
    python phase0b_nsp_klempt_connection.py
    python phase0b_nsp_klempt_connection.py --save
"""

import argparse
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

_HERE  = Path(__file__).resolve().parent
_NIFE  = Path("/home/nishioka/IKM_Hiwi/data_5species/main")
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_NIFE))

from hamilton_ode_jax_nsp import simulate_0d_nsp, theta_to_matrices

# ---------------------------------------------------------------------------
# Single-species NSP with nutrient modulation
# Species: So=0, An=1, Vd=2, Fn=3, Pg=4
# theta = [A_upper(15), b_diag(5)] = 20 params
# ---------------------------------------------------------------------------
N_SP    = 5

def _make_single_species_theta(a11=2.0, b1=1.5):
    """theta for single-species (So only): A=[[a11,0,...]], b=[b1,0,...].
    theta_full = [A_upper(15), b_diag(5)] — upper triangle row-major."""
    # A_upper: upper triangle of 5x5 symmetric A
    # A[i,i] = a11 for i=0, 0 otherwise
    A = np.zeros((N_SP, N_SP))
    A[0, 0] = a11
    A_upper = []
    for i in range(N_SP):
        for j in range(i, N_SP):
            A_upper.append(A[i, j])
    b = np.zeros(N_SP)
    b[0] = b1
    return np.array(A_upper + list(b), dtype=float)


def nsp_single_species_eq(c_local, n_steps=500, dt=1e-3):
    """Run NSP ODE (So only) at fixed c_local, return equilibrium φ_So."""
    theta = jnp.array(_make_single_species_theta())
    phi_init = jnp.zeros(N_SP).at[0].set(0.3)

    # c_hamilton = c_local * 100: scale to NSP internal units
    c_hamilton = max(0.01, float(c_local)) * 100.0
    traj = simulate_0d_nsp(
        theta, N_SP,
        n_steps=n_steps, dt=dt,
        phi_init=phi_init,
        c_const=c_hamilton, alpha_const=100.0,
    )
    # traj shape: (n_steps+1, N_SP) — phibar (mean-field φ)
    return float(traj[-1, 0])   # φ_So at equilibrium


def run(save=False):
    print(f"\n{'='*60}")
    print("Phase 0b: NSP (single-species So) ↔ Klempt PDE接続確認")
    print(f"{'='*60}\n")

    # Load Klempt PDE fields (Phase 0a output)
    npy_dir = _HERE
    phi_k = np.load(npy_dir / "klempt_phi_final.npy")
    c_k   = np.load(npy_dir / "klempt_c_final.npy")
    alpha_k = np.load(npy_dir / "klempt_alpha_final.npy")
    Nx, Ny = phi_k.shape

    print(f"  Klempt PDE fields: Nx={Nx} Ny={Ny}")
    print(f"  φ_total range : [{phi_k.min():.3f}, {phi_k.max():.3f}]")
    print(f"  c range       : [{c_k.min():.3f}, {c_k.max():.3f}]")
    print(f"  α range       : [{alpha_k.min():.3f}, {alpha_k.max():.3f}]")

    # Test NSP at several c_local values spanning Klempt range
    c_test_vals = [0.0, 0.01, 0.1, 0.5, 1.0]
    print(f"\n  NSP equilibrium φ_So vs c_local:")
    print(f"  {'c_local':>10}  {'φ_So_eq':>10}")
    phi_so_vals = []
    for c_loc in c_test_vals:
        phi_so = nsp_single_species_eq(c_loc)
        phi_so_vals.append(phi_so)
        print(f"  {c_loc:>10.3f}  {phi_so:>10.4f}")

    # Nutrient-dependence check: φ_So at c=1 vs c=0 should differ (even weakly)
    # NSP c_const modulates interaction strength — effect is subtle but monotonic
    nutrient_dependent = phi_so_vals[-1] > phi_so_vals[0]  # monotone increase
    all_positive = all(v >= 0 for v in phi_so_vals)
    equilibrium_positive = phi_so_vals[-1] > 0.0  # at c=1 (high nutrient)

    # Map NSP equilibrium over Klempt spatial field (downsampled 10×10 for speed)
    print(f"\n  Computing φ_So(x,y) map over Klempt c(x,y)...")
    step = max(1, Nx // 10)
    x_pts = np.arange(0, Nx, step)
    y_pts = np.arange(0, Ny, step)
    phi_so_map = np.zeros((len(x_pts), len(y_pts)))
    for i, xi in enumerate(x_pts):
        for j, yj in enumerate(y_pts):
            phi_so_map[i, j] = nsp_single_species_eq(float(c_k[xi, yj]))

    print(f"  φ_So map: [{phi_so_map.min():.4f}, {phi_so_map.max():.4f}]")
    # Where nutrient is low (biofilm interior), So should be suppressed
    c_low  = c_k[c_k < 0.1]
    c_high = c_k[c_k > 0.5]

    passed = equilibrium_positive and nutrient_dependent

    print(f"\n{'='*60}")
    print(f"BENCHMARK Phase 0b: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"  φ_So > 0 at high nutrient (c=1): {equilibrium_positive}  ({phi_so_vals[-1]:.4f})")
    print(f"  Nutrient-dependent dynamics     : {nutrient_dependent}")
    print(f"  φ_So map range over Klempt c(x,y): [{phi_so_map.min():.4f}, {phi_so_map.max():.4f}]")
    print(f"\n  Architectural note:")
    print(f"  NSP ODE governs WHICH species dominate (φᵢ fractions)")
    print(f"  Klempt PDE governs WHERE biofilm grows (φ_total density)")
    print(f"  Connection: φᵢ(x,y) = φ_total(x,y) × NSP_fraction_i(c(x,y))")
    print(f"{'='*60}\n")

    if save:
        _plot(phi_k, c_k, alpha_k, phi_so_map, x_pts, y_pts, Nx, Ny,
              c_test_vals, phi_so_vals, save)

    return passed, phi_so_map


def _plot(phi_k, c_k, alpha_k, phi_so_map, x_pts, y_pts, Nx, Ny,
          c_test_vals, phi_so_vals, save):
    import matplotlib.pyplot as plt
    import thesis_style
    fw, fh = thesis_style.use(width_frac=1.0, aspect=1.05)
    fig, axes = plt.subplots(2, 3, figsize=(fw, fh))

    # Klempt fields
    for ax, data, title, cmap in [
        (axes[0, 0], phi_k,   r"$\phi_\mathrm{total}(x,y)$",  "YlOrBr"),
        (axes[0, 1], c_k,     r"$c(x,y)$",                    "Blues"),
        (axes[0, 2], alpha_k, r"$\alpha(x,y)$",               "Greens"),
    ]:
        im = ax.imshow(data.T, origin="lower", cmap=cmap)
        ax.set_title(title)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

    # NSP φ_So(x,y) map over Klempt c field (downsampled)
    ax_map = axes[1, 0]
    im_map = ax_map.imshow(phi_so_map.T, origin="lower", cmap="Oranges",
                           extent=[0, Nx, 0, Ny])
    ax_map.set_title(r"$\phi_{So}(x,y)$ from NSP ($10\times10$ sample)")
    plt.colorbar(im_map, ax=ax_map, fraction=0.046, pad=0.02)

    # NSP φ_So vs c_local curve
    ax_curve = axes[1, 1]
    ax_curve.plot(c_test_vals, phi_so_vals, "o-", color="#d62728", lw=1.0, ms=3)
    ax_curve.set_xlabel(r"$c_\mathrm{local}$ (Klempt nutrient)")
    ax_curve.set_ylabel(r"$\phi_{So}$ equilibrium")
    ax_curve.set_title(r"\textit{S.\,oralis} equilibrium vs.\ nutrient")
    ax_curve.grid(True, alpha=0.3)

    # Connection summary (LaTeX-compatible text, no fontfamily override)
    ax_arch = axes[1, 2]
    ax_arch.axis("off")
    ax_arch.text(0.05, 0.95,
        r"\textbf{Multi-scale connection (Phase 0b)}" + "\n\n"
        r"Klempt PDE $\to$ $\phi$, $c$, $\alpha$" + "\n"
        r"$\downarrow$  $c(x,y)$ drives NSP" + "\n"
        r"NSP ODE $\to$ $\phi_{So}(x,y)$" + "\n"
        r"$\downarrow$  combine" + "\n"
        r"FEM: $E(\phi_i)$, $\varepsilon_g(\alpha)$, $\sigma_{vm}$",
        transform=ax_arch.transAxes,
        fontsize=7, va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))

    fig.suptitle(
        r"Phase 0b: Hamilton 5-species ODE $\leftrightarrow$ Klempt PDE connection check",
        fontsize=9)
    fig.tight_layout()

    out = _HERE / "phase0b_nsp_klempt_connection.png"
    out.unlink(missing_ok=True)   # 古い PNG を先に削除
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    run(save=args.save)
