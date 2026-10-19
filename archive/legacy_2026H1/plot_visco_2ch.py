"""
plot_visco_2ch.py  —  Visualise 2-ch Prony UMAT results for CH vs DH.
Run locally: python plot_visco_2ch.py --save visco_2ch_results.pdf
"""
import argparse, json, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Load ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
DATA_PATH  = SCRIPT_DIR / "visco_2ch_results.json"

K_RATIO   = 6.44 ** (1.0 / 2.68)   # ≈ 2.004
TAU1, TAU2 = 5.01, 80.72
E_INF, E1, E2 = 0.1977, 0.4553, 0.3470
T_PERIOD   = 1000.0

# Analytic relaxation modulus (normalised, t > 0)
def G_norm(t):
    return E_INF + E1 * np.exp(-t / TAU1) + E2 * np.exp(-t / TAU2)


def load():
    with open(DATA_PATH) as f:
        raw = json.load(f)
    # Standardise keys
    ch_key = [k for k in raw if "ch" in k.lower()][0]
    dh_key = [k for k in raw if "dh" in k.lower()][0]
    return raw[ch_key], raw[dh_key]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", default="visco_2ch_results.pdf")
    args = parser.parse_args()

    ch, dh = load()
    t_ch  = np.array(ch["t"])
    t_dh  = np.array(dh["t"])
    z_ch  = np.array(ch["z_mm"])
    z_dh  = np.array(dh["z_mm"])

    S11_ch   = np.array(ch["S11"])   * 1e3    # MPa → kPa
    S11_dh   = np.array(dh["S11"])   * 1e3
    SMIS_ch  = np.array(ch["SMISES"]) * 1e3
    SMIS_dh  = np.array(dh["SMISES"]) * 1e3

    # Index of top element (max compressive lateral stress)
    top = -1

    fig = plt.figure(figsize=(12, 9))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.35)

    # ── (A) Mises stress vs time at top element ───────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.plot(t_ch, SMIS_ch[:, top], "b-o", ms=4, label=f"CH (α_max={ch['z_mm'][top]:.3f}mm)")
    ax_a.plot(t_dh, SMIS_dh[:, top], "r--s", ms=4, label="DH (α_CH/k_ratio)")
    # Analytic long-time envelope
    t_fine = np.linspace(0, T_PERIOD, 500)
    ax_a.axhline(SMIS_ch[-1, top], color="b", lw=0.8, ls=":", alpha=0.5, label="CH long-time")
    ax_a.axhline(SMIS_dh[-1, top], color="r", lw=0.8, ls=":", alpha=0.5, label="DH long-time")
    ax_a.set_xlabel("Time [s]")
    ax_a.set_ylabel("von Mises stress [kPa]")
    ax_a.set_title("(A) Mises stress — top element", fontsize=9)
    ax_a.legend(fontsize=7)
    ax_a.set_xlim(0, T_PERIOD)

    # ── (B) Ratio CH/DH Mises at top ─────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    # only where both are non-zero
    mask = (SMIS_dh[:, top] > 1e-10)
    ratio = SMIS_ch[mask, top] / SMIS_dh[mask, top]
    ax_b.plot(t_ch[mask], ratio, "k-o", ms=4)
    ax_b.axhline(K_RATIO, color="gray", ls="--", label=f"k_ratio = {K_RATIO:.2f}")
    ax_b.set_xlabel("Time [s]")
    ax_b.set_ylabel("Mises_CH / Mises_DH")
    ax_b.set_title("(B) CH/DH stress ratio (→ k_ratio)", fontsize=9)
    ax_b.legend(fontsize=7)
    ax_b.set_xlim(0, T_PERIOD)

    # ── (C) S11 profile vs z at several times (CH) ───────────────────────────
    ax_c = fig.add_subplot(gs[1, 0])
    n_fr = len(t_ch)
    frames_to_plot = [0, n_fr//4, n_fr//2, 3*n_fr//4, n_fr-1]
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(frames_to_plot)))
    for k, fr in enumerate(frames_to_plot):
        ax_c.plot(S11_ch[fr, :], z_ch, color=cmap[k], lw=1.5,
                  label=f"t={t_ch[fr]:.0f}s")
    ax_c.set_xlabel("S11 [kPa]")
    ax_c.set_ylabel("z [mm]")
    ax_c.set_title("(C) CH lateral stress profile", fontsize=9)
    ax_c.legend(fontsize=6, loc="upper left")

    # ── (D) Normalised relaxation: Mises vs G_norm(t) * Mises_0 (CH) ─────────
    ax_d = fig.add_subplot(gs[1, 1])
    # At t=1000s the growth is maximum. Approximate "reference" as Mises_1000 / G_norm(0)
    mises_ch_top = SMIS_ch[:, top]
    # Plot G_norm * Mises_final vs t (expected if relaxation only)
    # This is a rough check — growth ramps so it's not a clean relaxation curve
    if mises_ch_top[-1] > 0:
        ax_d.plot(t_ch, mises_ch_top / mises_ch_top[-1], "b-o", ms=4, label="CH normalised Mises")
        mises_dh_top = SMIS_dh[:, top]
        ax_d.plot(t_dh, mises_dh_top / mises_dh_top[-1], "r--s", ms=4, label="DH normalised Mises")
        # Analytic: for pure relaxation after step load, G_norm(t)/G_norm(T_PERIOD)
        ax_d.plot(t_fine, G_norm(t_fine)/G_norm(T_PERIOD), "k:", lw=1.5,
                  label=f"G_norm(t)/G_norm({T_PERIOD:.0f}s)")
    ax_d.set_xlabel("Time [s]")
    ax_d.set_ylabel("Normalised stress")
    ax_d.set_title("(D) Relaxation shape (CH & DH vs Prony)", fontsize=9)
    ax_d.legend(fontsize=7)
    ax_d.set_xlim(0, T_PERIOD)

    fig.suptitle(
        "2-channel Prony UMAT — Abaqus/Standard biofilm viscoelasticity\n"
        f"G₀=2kPa · τ₁={TAU1:.1f}s · τ₂={TAU2:.1f}s · "
        f"e∞={E_INF:.3f} · e₁={E1:.3f} · k_ratio={K_RATIO:.2f}",
        fontsize=9
    )

    out = SCRIPT_DIR / args.save
    plt.savefig(out, bbox_inches="tight", dpi=150)
    print(f"Saved: {out}")

    # Print quantitative summary
    print(f"\nSummary (t=1000s, top element):")
    print(f"  CH: Mises={SMIS_ch[-1, top]:.4f} kPa  S11={S11_ch[-1, top]:.4f} kPa")
    print(f"  DH: Mises={SMIS_dh[-1, top]:.4f} kPa  S11={S11_dh[-1, top]:.4f} kPa")
    if SMIS_dh[-1, top] > 0:
        print(f"  Ratio CH/DH: {SMIS_ch[-1, top]/SMIS_dh[-1, top]:.3f}  (k_ratio={K_RATIO:.3f})")
    print(f"\nLong-time stress fraction (G_norm→e_inf={E_INF:.3f}):")
    print(f"  CH Mises(1000s) / Mises(max): {SMIS_ch[-1, top]/(SMIS_ch[:, top].max()+1e-20):.3f}")
    print(f"  DH Mises(1000s) / Mises(max): {SMIS_dh[-1, top]/(SMIS_dh[:, top].max()+1e-20):.3f}")


if __name__ == "__main__":
    main()
