"""
phase_field_at2_2d_sparse.py
============================
2-D AT2 phase-field fracture using scipy.sparse (structured FD grid).
Same physics as phase_field_at2_2d.py but no FEniCSx dependency.

PDE (damage field d, scalar):
  -G_c·ℓ·Δd + (G_c/ℓ + 2H)·d = 2H     Neumann on all boundaries
  H(x,z,t) = max_{s≤t} ½E·α(x,z,s)²    irreversibility

Discretisation: 5-point FD on structured Nx×Nz grid.
2-D Laplacian via Kronecker product of 1-D Neumann operators.
System matrix = G_c·ℓ·LAP (fixed) + diag(G_c/ℓ + 2H) (updated each step).
Solve with scipy.sparse.linalg.spsolve (direct, exact for moderate grids).
"""
from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.sparse import diags, eye, kron, lil_matrix
from scipy.sparse.linalg import spsolve

# ────────────────────────────────────────────────────────────────────────────
# Physical parameters (SI)
# ────────────────────────────────────────────────────────────────────────────
E_BIO  = 1.0e3     # Young's modulus [Pa]
G_C    = 1.0e-5    # Fracture energy [J/m²]
ELL    = 5.0e-6    # Phase-field length scale [m]
W      = 100.0e-6  # Domain width [m]
L      = 60.0e-6   # Film thickness [m]

K_RATIO  = 6.44 ** (1.0 / 2.68)   # k_eff^b(CH)/k_eff^b(DH) ≈ 2.00
ALPHA_C  = math.sqrt(G_C / (E_BIO * ELL))   # critical eigenstrain ≈ 0.045
N_STEPS  = 120
K_EFF_CH = ALPHA_C / N_STEPS * 1.6
K_EFF_DH = K_EFF_CH / K_RATIO


# ────────────────────────────────────────────────────────────────────────────
# Sparse Laplacian helpers
# ────────────────────────────────────────────────────────────────────────────
def _lap1d_neumann(n: int, h: float):
    """1-D Laplacian/h² with Neumann BCs via ghost nodes."""
    L = lil_matrix((n, n))
    for i in range(n):
        L[i, i] = -2.0
        if i > 0:
            L[i, i - 1] = 1.0
        if i < n - 1:
            L[i, i + 1] = 1.0
    # Neumann: ghost = neighbour → double the off-diagonal at boundaries
    L[0, 1]    += 1.0   # left:  d[-1]=d[1]  → coeff of d[1] becomes 2
    L[n-1, n-2] += 1.0  # right: d[n]=d[n-2] → coeff of d[n-2] becomes 2
    return L.tocsr() / h**2


def build_laplacian(Nx: int, Nz: int, hx: float, hz: float):
    """
    2-D Laplacian on Nx×Nz grid (node ordering: row-major, x outer, z inner).
    Returns sparse (Nx*Nz, Nx*Nz) matrix.
    """
    Lx = _lap1d_neumann(Nx, hx)
    Lz = _lap1d_neumann(Nz, hz)
    Ix = eye(Nx, format="csr")
    Iz = eye(Nz, format="csr")
    return kron(Ix, Lz, format="csr") + kron(Lx, Iz, format="csr")


# ────────────────────────────────────────────────────────────────────────────
# Growth profile
# ────────────────────────────────────────────────────────────────────────────
def make_k_field(k_base: float, profile: str,
                 X: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """
    k_eff^b at each grid node.  X, Z: 2-D arrays (Nx, Nz).
    """
    if profile == "uniform":
        phi = np.ones_like(Z)
    elif profile == "pg_substrate":
        # Pg concentrated near substrate z=0, characteristic depth L/4
        phi = np.exp(-Z / (L / 4.0))
    elif profile == "pg_patches":
        # CH-like lateral patches (hot spots) for visual interest
        phi = 1.0 + 0.5 * np.sin(2 * np.pi * X / W * 3)
        phi = np.clip(phi, 0.1, None)
    else:
        raise ValueError(profile)
    phi = phi / phi.mean()      # normalise: mean k_eff^b = k_base
    return (k_base * phi).ravel()   # flat, row-major


# ────────────────────────────────────────────────────────────────────────────
# Staggered AT2 solve
# ────────────────────────────────────────────────────────────────────────────
def run(k_base: float, profile: str,
        Nx: int, Nz: int, n_steps: int = N_STEPS,
        label: str = "") -> dict:

    hx = W / (Nx - 1)
    hz = L / (Nz - 1)
    x1d = np.linspace(0, W, Nx)
    z1d = np.linspace(0, L, Nz)
    X, Z = np.meshgrid(x1d, z1d, indexing="ij")   # (Nx, Nz)

    k_field = make_k_field(k_base, profile, X, Z)  # (Nx*Nz,)
    N = Nx * Nz

    LAP = build_laplacian(Nx, Nz, hx, hz)
    A_fixed = -G_C * ELL * LAP                     # diffusion part (fixed)

    d   = np.zeros(N)
    H   = np.zeros(N)
    results = []

    substrate_mask = (Z < hz * 2).ravel()   # nodes near z=0

    for step in range(1, n_steps + 1):
        t = step * (ALPHA_C * 1.6 / n_steps) / k_base

        alpha   = k_field * t
        psi_e   = 0.5 * E_BIO * alpha**2
        H       = np.maximum(H, psi_e)

        # System: (A_fixed + diag(G_c/ℓ + 2H)) d = 2H
        diag_vals = G_C / ELL + 2.0 * H
        A = A_fixed + diags(diag_vals, format="csr")
        rhs = 2.0 * H

        d_new = spsolve(A, rhs)
        d = np.maximum(d_new, d)            # irreversibility
        d = np.clip(d, 0.0, 1.0)

        d_max  = d.max()
        d_base = d[substrate_mask].max() if substrate_mask.any() else 0.0

        results.append({
            "step": step, "t": t,
            "alpha_mean": float(alpha.mean()),
            "d_max": d_max,
            "d_base": d_base,
            "d_2d": d.reshape(Nx, Nz).copy(),
        })

        if label and step % 30 == 0:
            print(f"  [{label}] step {step:3d}  ᾱ={alpha.mean():.3f}"
                  f"  d_max={d_max:.3f}  d_base={d_base:.3f}")

        if d_base > 0.95:
            if label:
                print(f"  [{label}] substrate fully damaged at step {step}")
            break

    return {"label": label, "results": results,
            "X": X, "Z": Z, "d_final": d.reshape(Nx, Nz)}


def t_crit(results: list, threshold: float = 0.5) -> float | None:
    for r in results:
        if r["d_max"] > threshold:
            return r["t"]
    return None


# ────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--nx", type=int, default=50, help="grid points in x")
    parser.add_argument("--nz", type=int, default=30, help="grid points in z")
    parser.add_argument("--n-steps", type=int, default=N_STEPS)
    parser.add_argument("--save", default="phase_field_at2_2d_sparse.pdf")
    args = parser.parse_args()

    Nx, Nz = args.nx, args.nz

    print("=" * 62)
    print("AT2 Phase-Field 2-D: Biofilm (scipy.sparse FD)")
    print("=" * 62)
    print(f"  E={E_BIO:.0f} Pa  G_c={G_C:.1e} J/m²  ℓ={ELL*1e6:.0f} µm")
    print(f"  Domain {W*1e6:.0f}×{L*1e6:.0f} µm  grid {Nx}×{Nz}={Nx*Nz} nodes")
    print(f"  α_c={ALPHA_C:.4f}  k_ratio={K_RATIO:.2f}")
    print()

    cases = [
        (K_EFF_CH, "uniform",      "CH"),
        (K_EFF_DH, "uniform",      "DH-unif"),
        (K_EFF_DH, "pg_substrate", "DH-Pg"),
    ]
    all_res = {}
    for k, prof, lbl in cases:
        print(f"Running {lbl:8s} ...")
        all_res[lbl] = run(k, prof, Nx, Nz, args.n_steps, label=lbl)
        print()

    # ── Summary ──────────────────────────────────────────────────────────
    print("-" * 62)
    tc_ch  = t_crit(all_res["CH"]["results"])
    tc_dhu = t_crit(all_res["DH-unif"]["results"])
    tc_dhp = t_crit(all_res["DH-Pg"]["results"])

    for lbl, out in all_res.items():
        tc = t_crit(out["results"])
        dm = out["results"][-1]["d_max"]
        tc_str = f"{tc:.1f} s" if tc else "not reached"
        print(f"  {lbl:10s}  t_crit={tc_str:>12s}  d_max_final={dm:.3f}")

    print()
    if tc_ch and tc_dhu:
        print(f"  Effect 1 (uniform):  t_crit(DH)/t_crit(CH) = "
              f"{tc_dhu/tc_ch:.2f}   theory = {K_RATIO:.2f}")
    if tc_ch and tc_dhp:
        print(f"  Effect 2 (Pg@sub):   t_crit(DH-Pg)/t_crit(CH) = "
              f"{tc_dhp/tc_ch:.2f}   (<1 → local interface damage)")

    if not args.plot:
        return

    # ── Plots ─────────────────────────────────────────────────────────────
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fracs  = [0.30, 0.65, 1.0]
    labels = ["Early", "Mid", "Final"]
    fig = plt.figure(figsize=(14, 9))
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.30)

    x_um = np.linspace(0, W * 1e6, list(all_res.values())[0]["X"].shape[0])
    z_um = np.linspace(0, L * 1e6, list(all_res.values())[0]["Z"].shape[1])

    for row, (lbl, out) in enumerate(all_res.items()):
        res = out["results"]

        for col, (frac, stage) in enumerate(zip(fracs, labels)):
            idx     = max(0, min(len(res)-1, int(len(res)*frac)))
            d_2d    = res[idx]["d_2d"]        # (Nx, Nz)
            alpha_m = res[idx]["alpha_mean"]

            ax = fig.add_subplot(gs[row, col])
            # pcolormesh(x1d, z1d, C): C must be (Nz, Nx) = d_2d.T
            im = ax.pcolormesh(x_um, z_um, d_2d.T,
                               vmin=0, vmax=1, cmap="hot_r",
                               shading="gouraud")
            ax.set_xlim(0, W*1e6);  ax.set_ylim(0, L*1e6)
            ax.set_xlabel("x [µm]", fontsize=7)
            ax.set_ylabel("z [µm]", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_title(f"{lbl} — {stage}\n"
                         f"ᾱ={alpha_m:.3f}  d_max={d_2d.max():.2f}",
                         fontsize=8)
            plt.colorbar(im, ax=ax, label="d", fraction=0.04, pad=0.02)

    tc_str = (f"{tc_dhu/tc_ch:.2f}" if (tc_ch and tc_dhu) else "N/A")
    fig.suptitle(
        "AT2 Phase-Field 2-D — Biofilm growth fracture: CH vs DH\n"
        f"E={E_BIO:.0f} Pa · G_c={G_C:.1e} J/m² · ℓ={ELL*1e6:.0f} µm · "
        f"k_ratio={K_RATIO:.2f} · t_crit ratio(uniform)={tc_str}",
        fontsize=10
    )
    plt.savefig(args.save, bbox_inches="tight", dpi=150)
    print(f"\nSaved: {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
