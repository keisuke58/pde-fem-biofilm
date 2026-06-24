"""
phase3_5species_stress.py
=========================
Phase 3: 5菌種 TMCMC posterior → DI(x,y) → E(x,y) → FEM応力

Pipeline:
  multiscale_coupling_2d.py の出力 (2d_fields_{condition}.npz)
    E_final(x,y)      ← DI→E マッピング済 (5菌種 TMCMC MAP)
    alpha_monod(x,y)  ← 局所成長活性 c/(k_M+c) × phi_total
    biofilm_mask(x,y) ← 卵形バイオフィルム領域
  → eps_growth = alpha_monod × biofilm_mask (Phase 0a α のプロキシ)
  → solve_2d_fem → σ_vm(x,y)
  → 4条件比較 + 条件間変位比

Klempt 2024 との対応 (Phase 3 拡張):
  単一菌種 E  → DI依存 E(x,y)  [5菌種 TMCMC]
  単一菌種 α  → alpha_monod(x,y) [Monodファクタ × biofilm]
  単一条件   → 4条件 (CS/CH/DS/DH) 比較

合格基準:
  条件間最大変位比 dysbiotic/commensal ≥ 2×
  (論文目標 > 5×, Phase 3はまず 2× を確認)

Run:
    # 前提: multiscale_coupling_2d.py が完了済 (FEM/_multiscale_2d_results/)
    python phase3_5species_stress.py
    python phase3_5species_stress.py --save
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import zoom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys; _sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import matplotlib.pyplot as plt
import thesis_style; fw, fh = thesis_style.use(width_frac=1.0, aspect=1.5)
import matplotlib.gridspec as gridspec

_HERE  = Path(__file__).resolve().parent
_FEM   = _HERE.parent
_MSCL  = _FEM / "_multiscale_2d_results"

sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_FEM))

from solve_stress_2d import solve_2d_fem

CONDITIONS = ["commensal_static", "commensal_hobic", "dysbiotic_static", "dysbiotic_hobic"]
LABELS     = ["Commensal\nStatic", "Commensal\nHOBiC", "Dysbiotic\nStatic", "Dysbiotic\nHOBiC"]
COLORS     = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]

E_MAX_PA  = 1000.0  # DI=0 (diverse/commensal)  [Pa]
E_MIN_PA  =   10.0  # DI=1 (single-species)       [Pa]
E_MATRX   =  900.0  # surrounding matrix (no biofilm) [Pa]


def di_to_E(di):
    """E = E_max*(1-DI)^2 + E_min*DI  (Shannon DI ∈ [0,1])."""
    return E_MAX_PA * (1 - di)**2 + E_MIN_PA * di


def load_condition(cond):
    """Load 2D multiscale fields for one condition."""
    npz = np.load(_MSCL / f"2d_fields_{cond}.npz")
    return {k: npz[k] for k in npz.keys()}


def load_0d_ref(cond):
    """Load 0D NSP equilibrium + DI for one condition."""
    import json
    rpath = _MSCL / f"ref_0d_{cond}.json"
    return json.load(open(rpath))


def resize_to(arr, target_shape):
    """Resize 2D array to target shape using zoom."""
    factors = (target_shape[0] / arr.shape[0], target_shape[1] / arr.shape[1])
    return zoom(arr, factors, order=1)


def run(save=False):
    print(f"\n{'='*60}")
    print("Phase 3: 5-species TMCMC → DI(x,y) → E(x,y) → FEM stress")
    print(f"{'='*60}\n")

    # Load Klempt α as eigenstrain source (Phase 0a output, 50×50)
    klempt_alpha = np.load(_HERE / "klempt_alpha_final.npy")

    results = {}
    for cond in CONDITIONS:
        data   = load_condition(cond)
        ref_0d = load_0d_ref(cond)
        Nx, Ny = data["biofilm_mask"].shape

        bio_mask = data["biofilm_mask"]      # (Nx, Ny) egg-shape
        di_0d    = float(ref_0d["di_0d"])    # scalar DI from 0D NSP posterior

        # E field: uniform DI_0d inside biofilm, matrix E outside
        E_bio    = di_to_E(di_0d)
        E_field  = np.where(bio_mask > 0.5, E_bio, E_MATRX)

        # ε_growth: resize Klempt α (50×50) to grid size (Nx×Ny)
        alpha_r  = resize_to(klempt_alpha, (Nx, Ny))
        eps_growth = alpha_r / 3.0 * bio_mask   # grow only inside biofilm

        print(f"  [{cond}]  DI_0d={di_0d:.3f}  E_bio={E_bio:.1f} Pa"
              f"  ε_g_max={eps_growth.max():.4f}")

        fem = solve_2d_fem(
            E_field          = E_field,
            nu               = 0.30,
            eps_growth_field = eps_growth,
            Nx               = Nx,
            Ny               = Ny,
            Lx               = 1.0,
            Ly               = 1.0,
            bc_type          = "bottom_fixed",
            stress_type      = "plane_strain",
        )

        u_mag = np.linalg.norm(fem["u"], axis=1).reshape(Nx, Ny)

        results[cond] = {
            "E_field":    E_field,
            "di_0d":      di_0d,
            "eps_growth": eps_growth,
            "bio_mask":   bio_mask,
            "sigma_vm":   fem["sigma_vm"].reshape(Nx-1, Ny-1),
            "u_mag":      u_mag,
            "u_max":      float(u_mag.max()),
            "s_max":      float(fem["sigma_vm"].max()),
        }

    # --- Condition comparison ---
    print(f"\n{'─'*60}")
    print(f"{'Condition':<25} {'σ_vm_max [Pa]':>14} {'|u|_max [m]':>12}")
    print(f"{'─'*60}")
    for cond, lbl in zip(CONDITIONS, LABELS):
        r = results[cond]
        print(f"  {cond:<23} {r['s_max']:>14.4f} {r['u_max']:>12.6f}")

    # Pass/Fail: dysbiotic / commensal displacement ratio ≥ 2
    u_cs = results["commensal_static"]["u_max"]
    u_ds = results["dysbiotic_static"]["u_max"]
    u_ch = results["commensal_hobic"]["u_max"]
    u_dh = results["dysbiotic_hobic"]["u_max"]
    ratio_static = u_ds / u_cs if u_cs > 1e-12 else float("inf")
    ratio_hobic  = u_dh / u_ch if u_ch > 1e-12 else float("inf")
    passed = (ratio_static >= 2.0) or (ratio_hobic >= 2.0)

    print(f"\n{'='*60}")
    print(f"BENCHMARK Phase 3: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"  |u| ratio  Dysbiotic/Commensal (Static) : {ratio_static:.2f}×")
    print(f"  |u| ratio  Dysbiotic/Commensal (HOBiC)  : {ratio_hobic:.2f}×")
    print(f"  Target: ≥2× (thesis goal: >5× with higher resolution)")
    print(f"{'='*60}\n")

    _plot(results, save)
    return results


def _plot(results, save):
    n = len(CONDITIONS)
    fig, axes = plt.subplots(3, n, figsize=(fw, fh))

    for col, (cond, lbl, color) in enumerate(zip(CONDITIONS, LABELS, COLORS)):
        r = results[cond]

        # Row 0: E field (DI-derived stiffness)
        im0 = axes[0, col].imshow(r["E_field"].T, origin="lower",
                                   cmap="RdYlGn", vmin=10, vmax=1000)
        axes[0, col].set_title(f"{lbl}\n(DI={r['di_0d']:.3f}, E_bio={di_to_E(r['di_0d']):.0f}Pa)",
                                fontsize=9, color=color, fontweight="bold")
        if col == 0:
            axes[0, col].set_ylabel(r"$E$(DI) [Pa]", fontsize=9)
        plt.colorbar(im0, ax=axes[0, col], fraction=0.046)

        # Row 1: σ_vm
        im1 = axes[1, col].imshow(r["sigma_vm"].T, origin="lower", cmap="hot")
        if col == 0:
            axes[1, col].set_ylabel(r"$\sigma_{vm}$ [Pa]", fontsize=9)
        plt.colorbar(im1, ax=axes[1, col], fraction=0.046)

        # Row 2: |u| displacement
        im2 = axes[2, col].imshow(r["u_mag"].T, origin="lower", cmap="cool")
        axes[2, col].set_xlabel("x")
        if col == 0:
            axes[2, col].set_ylabel(r"$|u|$ [m]", fontsize=9)
        plt.colorbar(im2, ax=axes[2, col], fraction=0.046)
        axes[2, col].set_title(
            f"|u|_max={r['u_max']:.5f} m\nσ_max={r['s_max']:.2f} Pa",
            fontsize=8,
        )

    fig.suptitle(
        "Phase 3: 5-species TMCMC posterior → DI(x,y) → E(x,y) → FEM stress\n"
        "Row 1: Young's modulus (DI-derived)  |  Row 2: von Mises stress  |  Row 3: displacement",
        fontsize=11,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if save:
        out = _HERE / "phase3_5species_stress.png"
        out.unlink(missing_ok=True)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    run(save=args.save)
