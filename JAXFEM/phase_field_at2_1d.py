"""
phase_field_at2_1d.py
=====================
1-D AT2 phase-field fracture for biofilm growth on a rigid substrate.

Physics
-------
Biofilm (thickness L, depth z∈[0,L]) grows laterally on a substrate.
Lateral constraint ⟹ biaxial eigenstress  σ₀ = −E·α(z,t)  (compression).
The stored elastic energy density

    Ψ_e(z, t) = ½ E α(z,t)²

acts as the driving force for damage d(z)∈[0,1] through the AT2 functional

    𝒻(d) = ∫₀ᴸ [ (1−d)² Ψ_e  +  G_c/(4ℓ) (d²/ℓ + ℓ (d')²) ] dz

Stationarity w.r.t. d gives the strong form (+ irreversibility via history H):

    −G_c·ℓ·d'' + (G_c/ℓ + 2H)·d = 2H,    d'(0)=d'(L)=0   [Neumann]

where  H(z,t) = max_{s≤t} Ψ_e(z,s)  enforces crack-irreversibility.

Growth model (from thesis)
--------------------------
    α(z,t) = k_eff^b(z) · t
    k_eff^b(z) = Σ_i φ_i(z) K_i^b        (composition-weighted growth coupling)

Link to FEM result:  σ̂ ∝ (k_eff^b)^2.68  →  ratio k(CH)/k(DH) = 6.44^(1/2.68) ≈ 2.00

Prediction of AT2 model
-----------------------
Critical eigenstrain for significant damage (d→0.5):  α_c = √(G_c / (E·ℓ))
Critical time:  t_crit = α_c / k_eff^b   →   t_crit(DH)/t_crit(CH) ≈ k(CH)/k(DH) ≈ 2

CH reaches significant damage ≈ 2× earlier than DH, consistent with σ_CH/σ_DH = 6.44×.

Parameters (calibrated to soft biofilm; all SI)
-----------------------------------------------
  E   = 1e3  Pa      Young's modulus (~1 kPa, AFM literature)
  G_c = 1e-5 J/m²   Critical energy release rate (~10 µJ/m², soft biofilm adhesion)
  ℓ   = 5e-6 m      Phase-field length scale (~1 cell diameter)
  L   = 60e-6 m     Film thickness

Usage
-----
  python phase_field_at2_1d.py               # print summary
  python phase_field_at2_1d.py --plot        # show/save figures
  python phase_field_at2_1d.py --n 200       # finer mesh
"""
from __future__ import annotations

import argparse
import sys
from typing import List, Tuple

import numpy as np
from scipy.linalg import solve_banded

# ---------------------------------------------------------------------------
# Physical parameters (SI)
# ---------------------------------------------------------------------------
E_BIO   = 1.0e3    # Young's modulus [Pa]
G_C     = 1.0e-5   # Fracture energy [J/m²]
ELL     = 5.0e-6   # Phase-field length scale [m]
L_BIO   = 60.0e-6  # Film thickness [m]

# Growth coupling ratio (from thesis power-law surrogate, MAP)
# k_eff^b(CH) / k_eff^b(DH) = 6.44^(1/2.68) ≈ 2.00
K_RATIO = 6.44 ** (1.0 / 2.68)   # ≈ 2.00

# Absolute growth rates: choose so CH damages within O(100) steps
# α_c = sqrt(G_c / (E·ℓ)) ≈ sqrt(1e-5 / (1e3 × 5e-6)) = sqrt(2) → set K so α=α_c in ~50 steps
ALPHA_C  = np.sqrt(G_C / (E_BIO * ELL))     # critical eigenstrain ≈ √2 ≈ 1.41
N_STEPS  = 120                               # growth increments
K_EFF_CH = ALPHA_C / N_STEPS * 1.5          # reach 1.5×α_c in N_STEPS steps
K_EFF_DH = K_EFF_CH / K_RATIO               # slower

# Depth profile shape (Pg concentration, from thesis FISH data)
# CH: uniform (Pg spread throughout depth)
# DH: Pg concentrated near substrate (z=0)
def phi_profile_ch(z: np.ndarray) -> np.ndarray:
    """Uniform Pg fraction → uniform k_eff^b."""
    return np.ones_like(z)  # relative, normalised to mean=1

def phi_profile_dh(z: np.ndarray) -> np.ndarray:
    """Pg concentrated near substrate (exponential drop-off)."""
    # More Pg near z=0 (substrate); characteristic depth ~ L/4
    return np.exp(-z / (L_BIO / 4.0))  # relative, peak at z=0

# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
def make_mesh(N: int, L: float) -> Tuple[np.ndarray, float]:
    z = np.linspace(0.0, L, N + 1)
    h = L / N
    return z, h

# ---------------------------------------------------------------------------
# Growth eigenstrain
# ---------------------------------------------------------------------------
def eigenstrain(z: np.ndarray, k_eff_base: float, phi_rel: np.ndarray,
                t: float) -> np.ndarray:
    """
    α(z, t) = k_eff^b(z) · t
    k_eff^b(z) = k_eff_base · phi_rel(z) / mean(phi_rel)
    (normalise so that mean k_eff^b = k_eff_base)
    """
    k_z = k_eff_base * phi_rel / phi_rel.mean()
    return k_z * t

# ---------------------------------------------------------------------------
# AT2 damage solve (1-D finite-difference, tridiagonal)
# ---------------------------------------------------------------------------
def solve_damage_fd(H: np.ndarray, d_prev: np.ndarray,
                    h: float, G_c: float, ell: float) -> np.ndarray:
    """
    Solve  -G_c·ℓ·d'' + (G_c/ℓ + 2H)·d = 2H
    with Neumann BCs d'(0) = d'(L) = 0, using second-order FD.

    Discretisation (central differences, ghost nodes for Neumann):
      Interior i∈[1,N-1]:
        -c_d·(d[i+1] - 2d[i] + d[i-1])/h² + r[i]·d[i] = 2H[i]
      Boundary i=0: d[-1]=d[1]  →  -2c_d·(d[1]-d[0])/h² + r[0]·d[0] = 2H[0]
      Boundary i=N: d[N+1]=d[N-1] (same)

    Returns d clipped to [d_prev, 1] (irreversibility).
    """
    n = len(H)
    c_d = G_c * ell          # diffusion coefficient
    r   = G_c / ell + 2.0 * H  # diagonal reaction term

    # Assemble tridiagonal: ab[0]=superdiag, ab[1]=main, ab[2]=subdiag
    ab = np.zeros((3, n))
    rhs = 2.0 * H.copy()

    # Interior nodes
    for i in range(1, n - 1):
        ab[1, i] = c_d * 2.0 / h**2 + r[i]
        ab[0, i] = -c_d / h**2          # superdiag at (i, i+1)
        ab[2, i] = -c_d / h**2          # subdiag at (i, i-1)

    # Left boundary i=0: ghost node d[-1] = d[1]  → coefficient of d[1] doubles
    ab[1, 0] = c_d * 2.0 / h**2 + r[0]
    ab[0, 0] = -c_d * 2.0 / h**2       # doubled

    # Right boundary i=n-1: ghost node d[n] = d[n-2] → coefficient of d[n-2] doubles
    ab[1, n-1] = c_d * 2.0 / h**2 + r[n-1]
    ab[2, n-1] = -c_d * 2.0 / h**2     # doubled

    # Fix superdiag offset: scipy solve_banded stores ab[0, j] = A[i-1, i]
    # i.e. the superdiag entry for column j is ab[0, j].
    # Our ab[0, i] is A[i, i+1] → shift by 1 for scipy convention.
    ab_scipy = np.zeros((3, n))
    ab_scipy[1, :] = ab[1, :]          # main diagonal
    ab_scipy[0, 1:] = ab[0, :-1]       # superdiag: ab[0,i] → ab_scipy[0, i+1]
    ab_scipy[2, :-1] = ab[2, 1:]       # subdiag:   ab[2,i] → ab_scipy[2, i-1]

    d = solve_banded((1, 1), ab_scipy, rhs)
    d = np.clip(d, d_prev, 1.0)
    return d

# ---------------------------------------------------------------------------
# Main simulation loop
# ---------------------------------------------------------------------------
def run(k_eff: float, phi_rel: np.ndarray, z: np.ndarray, h: float,
        n_steps: int = N_STEPS, label: str = "") -> List[dict]:
    """
    Staggered AT2 growth simulation.
    Returns list of result dicts per step.
    """
    n = len(z)
    d = np.zeros(n)
    H = np.zeros(n)
    results = []

    for step in range(1, n_steps + 1):
        t = step * (ALPHA_C * 1.5 / n_steps) / k_eff   # physical time [s]
        alpha = eigenstrain(z, k_eff, phi_rel, t)
        psi_e = 0.5 * E_BIO * alpha**2
        H = np.maximum(H, psi_e)
        d = solve_damage_fd(H, d, h, G_C, ELL)

        results.append({
            "step": step,
            "t": t,
            "alpha_mean": alpha.mean(),
            "d_max": d.max(),
            "d_base": d[0],          # damage at substrate (z=0)
            "sigma_mean": -(1 - d)**2 * E_BIO * alpha,   # compressive (negative)
            "d": d.copy(),
            "alpha": alpha.copy(),
        })

        if d[0] > 0.95:
            if label:
                print(f"  [{label}] substrate fully damaged at step {step}")
            break

    return results


def t_crit(results: List[dict], threshold: float = 0.5) -> float | None:
    """Return time when d_max first exceeds threshold."""
    for r in results:
        if r["d_max"] > threshold:
            return r["t"]
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot", action="store_true", help="Show/save figures")
    parser.add_argument("-n", "--n-elem", type=int, default=120,
                        help="Number of FD elements (default 120)")
    parser.add_argument("--n-steps", type=int, default=N_STEPS,
                        help="Growth increments (default %(default)s)")
    parser.add_argument("--save", type=str, default="phase_field_at2_1d.pdf",
                        help="Output PDF name")
    args = parser.parse_args()

    N = args.n_elem
    z, h = make_mesh(N, L_BIO)
    z_um = z * 1e6  # µm for plotting

    phi_ch = phi_profile_ch(z)
    phi_dh = phi_profile_dh(z)

    print("=" * 60)
    print("AT2 Phase-Field 1-D: Biofilm Growth Fracture")
    print("=" * 60)
    print(f"  E    = {E_BIO:.0f} Pa")
    print(f"  G_c  = {G_C:.2e} J/m²")
    print(f"  ℓ    = {ELL*1e6:.1f} µm")
    print(f"  L    = {L_BIO*1e6:.0f} µm")
    print(f"  α_c  = {ALPHA_C:.3f}  (critical eigenstrain, d→0.5)")
    print(f"  k_eff^b ratio (CH/DH): {K_RATIO:.2f}")
    print()

    print(f"Running CH      (k_eff^b = {K_EFF_CH:.4e}, uniform depth) ...")
    res_ch = run(K_EFF_CH, phi_ch, z, h, n_steps=args.n_steps, label="CH")

    phi_dh_unif = phi_profile_ch(z)   # same shape (uniform) but lower k_eff
    print(f"Running DH-unif (k_eff^b = {K_EFF_DH:.4e}, uniform depth) ...")
    res_dh_u = run(K_EFF_DH, phi_dh_unif, z, h, n_steps=args.n_steps, label="DH-unif")

    print(f"Running DH-Pg   (k_eff^b = {K_EFF_DH:.4e}, Pg@substrate) ...")
    res_dh = run(K_EFF_DH, phi_dh, z, h, n_steps=args.n_steps, label="DH-Pg")

    tc_ch   = t_crit(res_ch)
    tc_dhu  = t_crit(res_dh_u)
    tc_dh   = t_crit(res_dh)
    ratio_unif = f"{tc_dhu/tc_ch:.2f}" if (tc_ch and tc_dhu) else "N/A"
    ratio_pg   = f"{tc_dh/tc_ch:.2f}"  if (tc_ch and tc_dh)  else "N/A"
    print()
    print("  Effect 1 — uniform profiles (k_eff^b controls):")
    print(f"    t_crit(CH)      = {tc_ch:.1f} s" if tc_ch else "    t_crit(CH) = not reached")
    print(f"    t_crit(DH-unif) = {tc_dhu:.1f} s" if tc_dhu else "    t_crit(DH-unif) = not reached")
    print(f"    t_crit(DH-unif)/t_crit(CH) = {ratio_unif}  (theory = k_ratio = {K_RATIO:.2f})")
    print()
    print("  Effect 2 — DH with Pg depth concentration:")
    print(f"    t_crit(DH-Pg)   = {tc_dh:.1f} s" if tc_dh else "    t_crit(DH-Pg) = not reached")
    print(f"    t_crit(DH-Pg)/t_crit(CH) = {ratio_pg}  (< 1 → local substrate damage faster)")
    print()
    print("  Physical interpretation:")
    print("    CH: uniform eigenstress → bulk damage (σ_CH > σ_DH globally)")
    print("    DH: Pg@substrate amplifies local k_eff^b → interfacial damage concentrates")
    print("        → 2 distinct failure modes despite lower global stress in DH")

    if not args.plot:
        return

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # Three snapshots: ~30%, ~60%, final
    def snap_idx(results, frac):
        return max(0, min(len(results) - 1, int(len(results) * frac)))

    fracs = [0.30, 0.65, 1.0]
    titles = ["Early growth", "Intermediate", "Final"]
    col_ch = "#2171B5"
    col_dh = "#CB181D"

    for col, (frac, title) in enumerate(zip(fracs, titles)):
        idx_ch = snap_idx(res_ch, frac)
        idx_dh = snap_idx(res_dh, frac)
        r_ch = res_ch[idx_ch]
        r_dh = res_dh[idx_dh]

        # Damage profile
        ax = fig.add_subplot(gs[0, col])
        ax.plot(z_um, r_ch["d"], color=col_ch, lw=2.0, label="CH")
        ax.plot(z_um, r_dh["d"], color=col_dh, lw=2.0, ls="--", label="DH")
        ax.axvline(0, color="k", lw=0.8, ls=":")
        ax.set_xlim(0, L_BIO * 1e6)
        ax.set_ylim(-0.05, 1.10)
        ax.set_xlabel("depth z [µm]", fontsize=9)
        ax.set_ylabel("damage d", fontsize=9)
        ax.set_title(f"{title}\n"
                     f"CH: ᾱ={r_ch['alpha_mean']:.3f}  "
                     f"DH: ᾱ={r_dh['alpha_mean']:.3f}", fontsize=8)
        ax.legend(fontsize=8, loc="upper right")

        # Stress profile (compressive → negative; show magnitude)
        ax2 = fig.add_subplot(gs[1, col])
        sig_ch = np.abs(r_ch["sigma_mean"]) * 1e-3   # kPa
        sig_dh = np.abs(r_dh["sigma_mean"]) * 1e-3
        ax2.plot(z_um, sig_ch, color=col_ch, lw=2.0, label="CH")
        ax2.plot(z_um, sig_dh, color=col_dh, lw=2.0, ls="--", label="DH")
        ax2.set_xlim(0, L_BIO * 1e6)
        ax2.set_xlabel("depth z [µm]", fontsize=9)
        ax2.set_ylabel("|σ| [kPa]", fontsize=9)
        ax2.legend(fontsize=8)

    fig.suptitle(
        "AT2 Phase-Field (1-D) — Biofilm Growth Fracture: CH vs DH\n"
        f"E={E_BIO:.0f} Pa, G_c={G_C:.1e} J/m², ℓ={ELL*1e6:.0f} µm  |  "
        f"k_ratio(CH/DH)={K_RATIO:.2f}  →  "
        r"$t_\mathrm{crit}$(DH)/$t_\mathrm{crit}$(CH) ≈ " + ratio_str,
        fontsize=11
    )

    plt.savefig(args.save, bbox_inches="tight", dpi=150)
    print(f"Saved {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
