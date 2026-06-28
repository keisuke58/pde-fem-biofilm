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
    Profiles:
      z_gradient   - growth ∝ (z/L)^0.8 from 1D Monod model (max at free surface)
      z_patchy     - z_gradient × lateral sine patches (CH heterogeneous colonisation)
      pg_substrate - Pg-only, concentrated near z=0
      pg_combo     - z_gradient (DH baseline) + Pg substrate component
      uniform      - spatially constant (kept for reference)
    """
    z_norm = Z / L   # 0=substrate, 1=free surface
    z_profile = z_norm ** 0.8 + 0.01   # from 1D bar alpha profile; +0.01 avoids exact zero

    if profile == "z_gradient":
        # commensal-dominated: aerobic bacteria grow best near free surface
        phi = z_profile
    elif profile == "z_patchy":
        # CH patchy colonisation: z-gradient × lateral heterogeneity
        patch = 1.0 + 0.6 * np.sin(2 * np.pi * X / W * 2.5)
        phi = z_profile * np.clip(patch, 0.1, None)
    elif profile == "pg_substrate":
        # Pg-only: concentrated near substrate
        phi = np.exp(-Z / (L / 4.0))
    elif profile == "pg_combo":
        # DH: reduced baseline z-gradient (1/k_ratio) + localised Pg at substrate
        phi_base = z_profile / K_RATIO
        phi_pg   = 0.8 * np.exp(-Z / (L / 5.0))
        phi = phi_base + phi_pg
    elif profile == "uniform":
        phi = np.ones_like(Z)
    else:
        raise ValueError(profile)
    phi = phi / phi.mean()   # normalise to mean=1 so k_base sets absolute rate
    return (k_base * phi).ravel()


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
        (K_EFF_CH, "z_gradient",   "CH"),      # commensal: z-gradient, fast (k_ratio×DH)
        (K_EFF_DH, "z_gradient",   "DH-base"), # dysbiotic baseline: same shape, 2× slower
        (K_EFF_DH, "pg_substrate", "DH-Pg"),   # dysbiotic + Pg: substrate damage, fastest
    ]
    all_res = {}
    for k, prof, lbl in cases:
        print(f"Running {lbl:8s} ...")
        all_res[lbl] = run(k, prof, Nx, Nz, args.n_steps, label=lbl)
        print()

    # ── Summary ──────────────────────────────────────────────────────────
    print("-" * 62)
    tc_ch  = t_crit(all_res["CH"]["results"])
    tc_dhb = t_crit(all_res["DH-base"]["results"])
    tc_dhp = t_crit(all_res["DH-Pg"]["results"])

    for lbl, out in all_res.items():
        tc = t_crit(out["results"])
        dm = out["results"][-1]["d_max"]
        tc_str = f"{tc:.1f} s" if tc else "not reached"
        print(f"  {lbl:10s}  t_crit={tc_str:>12s}  d_max_final={dm:.3f}")

    print()
    if tc_ch and tc_dhb:
        print(f"  Effect 1 (baseline): t_crit(DH)/t_crit(CH) = "
              f"{tc_dhb/tc_ch:.2f}   theory ≈ {K_RATIO:.2f}")
    if tc_ch and tc_dhp:
        print(f"  Effect 2 (DH+Pg):    t_crit(DH-Pg)/t_crit(CH) = "
              f"{tc_dhp/tc_ch:.2f}   (<1 → interface damage precedes bulk)")

    if not args.plot:
        return

    # ── Plots ─────────────────────────────────────────────────────────────
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    tc_dhb = t_crit(all_res["DH-base"]["results"])
    tc_dhp = t_crit(all_res["DH-Pg"]["results"])
    ratio_str = (f"{tc_dhb/tc_ch:.2f}" if (tc_ch and tc_dhb) else "N/A")
    pg_str    = (f"{tc_dhp/tc_ch:.2f}" if (tc_ch and tc_dhp) else "N/A")

    fig = plt.figure(figsize=(14, 10))
    # 3 rows: Row 0=CH heatmaps, Row 1=DH-Pg heatmaps, Row 2=time curves
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.55, wspace=0.32,
                           height_ratios=[1, 1, 0.8])

    x_um = np.linspace(0, W * 1e6, list(all_res.values())[0]["X"].shape[0])
    z_um = np.linspace(0, L * 1e6, list(all_res.values())[0]["Z"].shape[1])

    # ── Rows 0–1: 2D heatmaps for CH and DH-Pg (most distinct spatial patterns) ──
    map_cases = [("CH", 0), ("DH-Pg", 1)]
    fracs  = [0.30, 0.65, 1.0]
    stages = ["Early", "Mid", "Final"]

    for lbl, row in map_cases:
        out = all_res[lbl]
        res = out["results"]
        for col, (frac, stage) in enumerate(zip(fracs, stages)):
            idx    = max(0, min(len(res)-1, int(len(res)*frac)))
            d_2d   = res[idx]["d_2d"]      # (Nx, Nz)
            t_real = res[idx]["t"]
            dmax   = d_2d.max()
            ax = fig.add_subplot(gs[row, col])
            im = ax.pcolormesh(x_um, z_um, d_2d.T,
                               vmin=0, vmax=1, cmap="hot_r", shading="gouraud")
            ax.set_xlim(0, W*1e6);  ax.set_ylim(0, L*1e6)
            ax.set_xlabel("x [µm]", fontsize=7)
            ax.set_ylabel("z [µm]", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_title(f"{lbl} — {stage}  (t={t_real:.0f}s)\n"
                         f"d_max={dmax:.2f}", fontsize=8)
            plt.colorbar(im, ax=ax, label="d", fraction=0.04, pad=0.02)

    # ── Row 2: d_max(t) and d_base(t) time comparison ───────────────────────
    ax_l = fig.add_subplot(gs[2, :2])   # spans columns 0–1
    ax_r = fig.add_subplot(gs[2, 2])    # column 2

    colors  = {"CH": "blue", "DH-base": "green", "DH-Pg": "red"}
    markers = {"CH": "o",    "DH-base": "s",     "DH-Pg": "^"}
    for lbl, out in all_res.items():
        res = out["results"]
        ts    = [r["t"]      for r in res]
        dmaxs = [r["d_max"]  for r in res]
        dbas  = [r["d_base"] for r in res]
        c, m  = colors[lbl], markers[lbl]
        ax_l.plot(ts, dmaxs, color=c, marker=m, ms=3, lw=1.5, label=f"{lbl} d_max")
        ax_l.plot(ts, dbas,  color=c, marker=m, ms=3, lw=0.8, ls="--",
                  label=f"{lbl} d_base")

    ax_l.axhline(0.5, color="gray", ls=":", lw=0.8, label="d=0.5 (t_crit)")
    ax_l.set_xlabel("Time [s]")
    ax_l.set_ylabel("Damage d")
    ax_l.set_title("d_max(t) and d_substrate(t) — all cases", fontsize=8)
    ax_l.legend(fontsize=6, ncol=2, loc="upper left")
    ax_l.set_ylim(0, 1.05)

    # d_max bar chart at t_crit of CH
    bars = {"CH": all_res["CH"]["results"][-1]["d_max"],
            "DH-base": all_res["DH-base"]["results"][-1]["d_max"],
            "DH-Pg": all_res["DH-Pg"]["results"][-1]["d_max"]}
    ax_r.bar(list(bars.keys()), list(bars.values()),
             color=[colors[k] for k in bars], alpha=0.75)
    ax_r.set_ylabel("d_max at end of simulation")
    ax_r.set_title("Final damage state", fontsize=8)
    ax_r.set_ylim(0, 1.1)
    for i, (k, v) in enumerate(bars.items()):
        ax_r.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)

    fig.suptitle(
        "AT2 Phase-Field 2-D — Biofilm growth fracture: CH vs DH\n"
        f"CH (z-gradient, k_CH): surface damage, fast  |  "
        f"DH-Pg (Pg@substrate, k_DH): interface damage\n"
        f"E={E_BIO:.0f} Pa · G_c={G_C:.1e} J/m^2 · l={ELL*1e6:.0f} um · "
        f"k_ratio={K_RATIO:.2f} · t_crit(DH-base/CH)={ratio_str} · "
        f"t_crit(DH-Pg/CH)={pg_str}",
        fontsize=9
    )
    plt.savefig(args.save, bbox_inches="tight", dpi=150)
    print(f"\nSaved: {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
