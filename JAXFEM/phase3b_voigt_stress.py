"""
phase3b_voigt_stress.py
=======================
Phase 3b: Voigt則による種別E場 → FEM応力（Phase 3の生物学的改良版）

Phase 3 の問題点:
  DI基準 E(DI) = E_max*(1-DI)^2 + E_min*DI は
  CH(DI=0.843, So=94%→硬い) に E=33 Pa を割り当てる（生物学的逆転）。

Voigt則（加重平均）:
  E(x,y) = Σᵢ φᵢ(x,y) × E_sp_i
  E_SPECIES = [So=1000, An=800, Vd=600, Fn=200, Pg=10]  [Pa]
  So (S.oralis): 主要EPS産生菌。グルカン骨格＝剛性担当。
  Vd (Veillonella): 乳酸代謝→So競合。EPS産生ほぼゼロ→骨格崩壊。

生物学的ストーリー:
  CS/CH/DS: So 70-94% → EPS scaffold豊富 → 硬い (E≈800-960 Pa)
  DH:       Vd 47%    → EPS scaffold崩壊 → 軟らかい (E≈500 Pa)
  → 同じ固有ひずみα でも DH は大きな変位 → インプラント界面の機械的脆弱化

合格基準（Phase 3b）:
  Voigt則で |u|_DH / |u|_CH ≥ 1.5× （生物学的に正しい方向）
  且つ DH が全4条件中最大変位

Run:
    python phase3b_voigt_stress.py
    python phase3b_voigt_stress.py --save
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

_HERE = Path(__file__).resolve().parent
_FEM  = _HERE.parent
_MSCL = _FEM / "_multiscale_2d_results"

sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_FEM))

from solve_stress_2d import solve_2d_fem

CONDITIONS = ["commensal_static", "commensal_hobic", "dysbiotic_static", "dysbiotic_hobic"]
LABELS     = ["Commensal\nStatic", "Commensal\nHOBiC", "Dysbiotic\nStatic", "Dysbiotic\nHOBiC"]
COLORS     = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
SPECIES    = ["So", "An", "Vd", "Fn", "Pg"]

# Species-specific E [Pa] — EPS産生能 = 剛性代理指標
# So: グルカンEPS主産生菌 (Streptococcus oralis)
# An: 中程度EPS (Actinomyces naeslundii)
# Vd: EPS産生ほぼゼロ (Veillonella dispar) ← DHで47%
# Fn: 線維状付着のみ (Fusobacterium nucleatum)
# Pg: 病原性高いが剛性貢献低 (Porphyromonas gingivalis)
E_SPECIES_PA = np.array([1000.0, 800.0, 600.0, 200.0, 10.0])
E_MATRIX_PA  = 900.0  # 周囲基質 (インプラント周囲骨/軟組織)


def voigt_E(phi_vec):
    """Voigt則: E = Σ φᵢ × Eᵢ  (加重平均 = EPS産生能の線形混合)"""
    return float(np.dot(phi_vec, E_SPECIES_PA))


def resize_to(arr, target_shape):
    factors = (target_shape[0] / arr.shape[0], target_shape[1] / arr.shape[1])
    return zoom(arr, factors, order=1)


def run(save=False):
    print(f"\n{'='*65}")
    print("Phase 3b: Voigt則 E = Σφᵢ×Eᵢ → FEM応力")
    print(f"{'='*65}\n")

    klempt_alpha = np.load(_HERE / "klempt_alpha_final.npy")

    # --- 各条件のφᵢとE_Voigt表示 ---
    import json
    print(f"{'Condition':<22} {'DI':>5} | {'So':>5} {'An':>5} {'Vd':>5} {'Fn':>5} {'Pg':>5} |"
          f" {'E_DI':>7} {'E_Voigt':>8}")
    print("-" * 80)
    for cond in CONDITIONS:
        d = json.load(open(_MSCL / f"ref_0d_{cond}.json"))
        phi = np.array(d["phi_final"])
        di  = float(d["di_0d"])
        E_di    = 1000*(1-di)**2 + 10*di
        E_vgt   = voigt_E(phi)
        pcts    = "  ".join(f"{p*100:5.1f}" for p in phi)
        print(f"  {cond:<20} {di:>5.3f} | {pcts} | {E_di:>7.1f} {E_vgt:>8.1f}")
    print()

    results = {}
    for cond in CONDITIONS:
        d = json.load(open(_MSCL / f"ref_0d_{cond}.json"))
        phi_vec  = np.array(d["phi_final"])   # (5,)
        E_vgt    = voigt_E(phi_vec)

        data     = np.load(_MSCL / f"2d_fields_{cond}.npz")
        bio_mask = data["biofilm_mask"]
        Nx, Ny   = bio_mask.shape

        # E場: Voigt値をバイオフィルム内に均一配置、外は基質E
        E_field  = np.where(bio_mask > 0.5, E_vgt, E_MATRIX_PA)

        # 固有ひずみ: Klempt α を Voigt版でも使用
        alpha_r    = resize_to(klempt_alpha, (Nx, Ny))
        eps_growth = alpha_r / 3.0 * bio_mask

        print(f"  [{cond}]  E_bio={E_vgt:.1f} Pa  E_mat={E_MATRIX_PA:.0f} Pa  "
              f"ε_g_max={eps_growth.max():.4f}")

        fem = solve_2d_fem(
            E_field          = E_field,
            nu               = 0.30,
            eps_growth_field = eps_growth,
            Nx = Nx, Ny = Ny, Lx = 1.0, Ly = 1.0,
            bc_type          = "bottom_fixed",
            stress_type      = "plane_strain",
        )

        u_mag = np.linalg.norm(fem["u"], axis=1).reshape(Nx, Ny)
        results[cond] = {
            "phi_vec":   phi_vec,
            "E_voigt":   E_vgt,
            "E_field":   E_field,
            "sigma_vm":  fem["sigma_vm"].reshape(Nx-1, Ny-1),
            "u_mag":     u_mag,
            "u_max":     float(u_mag.max()),
            "s_max":     float(fem["sigma_vm"].max()),
        }

    # --- 結果比較 ---
    print(f"\n{'─'*65}")
    print(f"{'Condition':<22} {'E_bio[Pa]':>10} {'σ_max[Pa]':>10} {'|u|_max[m]':>12}")
    print(f"{'─'*65}")
    for cond, lbl in zip(CONDITIONS, LABELS):
        r = results[cond]
        print(f"  {cond:<20} {r['E_voigt']:>10.1f} {r['s_max']:>10.4f} {r['u_max']:>12.6f}")

    # Pass/Fail: DH最大変位 & DH/CH比 ≥ 1.5×
    u_ch = results["commensal_hobic"]["u_max"]
    u_dh = results["dysbiotic_hobic"]["u_max"]
    u_cs = results["commensal_static"]["u_max"]
    u_ds = results["dysbiotic_static"]["u_max"]
    ratio_hobic  = u_dh / u_ch if u_ch > 1e-12 else float("inf")
    ratio_static = u_ds / u_cs if u_cs > 1e-12 else float("inf")
    dh_is_max    = u_dh >= max(results[c]["u_max"] for c in CONDITIONS)
    passed       = (ratio_hobic >= 1.5) and dh_is_max

    print(f"\n  生物学的期待: DH(Vd支配・軟) > CH(So支配・硬) の変位")
    print(f"  |u|_DH / |u|_CH  = {ratio_hobic:.3f}×  (閾値 ≥1.5×)")
    print(f"  |u|_DS / |u|_CS  = {ratio_static:.3f}×")
    print(f"  DH が全4条件中最大変位: {dh_is_max}")

    print(f"\n{'='*65}")
    print(f"BENCHMARK Phase 3b: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"  Voigt則: E_CH={results['commensal_hobic']['E_voigt']:.0f}Pa"
          f" > E_DH={results['dysbiotic_hobic']['E_voigt']:.0f}Pa")
    print(f"  → 軟らかいDH biofilm = 同じαで大きな変位 = 機械的脆弱化")
    print(f"{'='*65}\n")

    if save:
        _plot(results, save)
    return results


def _plot(results, save):
    n = len(CONDITIONS)
    fig, axes = plt.subplots(3, n, figsize=(fw, fh))

    phi_all = np.array([results[c]["phi_vec"] for c in CONDITIONS])  # (4,5)

    for col, (cond, lbl, color) in enumerate(zip(CONDITIONS, LABELS, COLORS)):
        r = results[cond]

        # Row 0: 種組成 bar
        ax = axes[0, col]
        bars = ax.bar(SPECIES, r["phi_vec"]*100,
                      color=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"])
        ax.set_ylim(0, 100)
        ax.set_title(f"{lbl}\nE_bio={r['E_voigt']:.0f} Pa (Voigt)",
                     fontsize=9, color=color, fontweight="bold")
        if col == 0:
            ax.set_ylabel(r"$\phi_i$ [\%]", fontsize=9)
        ax.axhline(50, color="gray", lw=0.5, ls="--")

        # Row 1: σ_vm
        im1 = axes[1, col].imshow(r["sigma_vm"].T, origin="lower", cmap="hot")
        if col == 0:
            axes[1, col].set_ylabel(r"$\sigma_{vm}$ [Pa]", fontsize=9)
        plt.colorbar(im1, ax=axes[1, col], fraction=0.046)
        axes[1, col].set_title(r"$\sigma_{max}$=" + f"{r['s_max']:.2f} Pa", fontsize=8)

        # Row 2: |u|
        im2 = axes[2, col].imshow(r["u_mag"].T, origin="lower", cmap="cool")
        axes[2, col].set_xlabel("x")
        if col == 0:
            axes[2, col].set_ylabel(r"$|u|$ [m]", fontsize=9)
        plt.colorbar(im2, ax=axes[2, col], fraction=0.046)
        axes[2, col].set_title(r"$|u|_{max}$=" + f"{r['u_max']:.5f} m", fontsize=8)

    # 注釈: E_SPECIES_PA
    fig.text(0.01, 0.01,
             f"E_species [Pa]: So={E_SPECIES_PA[0]:.0f}  An={E_SPECIES_PA[1]:.0f}"
             f"  Vd={E_SPECIES_PA[2]:.0f}  Fn={E_SPECIES_PA[3]:.0f}  Pg={E_SPECIES_PA[4]:.0f}"
             f"  |  E_matrix={E_MATRIX_PA:.0f}\n"
             r"Voigt: $E_{bio} = \sum_i \phi_i \cdot E_i$ (EPS-weighted average)",
             fontsize=8, va="bottom", ha="left",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    fig.suptitle(
        "Phase 3b: Voigt-rule $E(\\phi_i)$ $\\rightarrow$ FEM stress\n"
        "So-dominated (CH/DS: stiff) vs Vd-dominated (DH: soft) -- Voigt homogenisation",
        fontsize=11,
    )
    plt.tight_layout(rect=[0, 0.07, 1, 0.95])

    if save:
        out = _HERE / "phase3b_voigt_stress.png"
        out.unlink(missing_ok=True)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    run(save=args.save)
