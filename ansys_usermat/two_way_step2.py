#!/usr/bin/env python3
"""Two-way coupling, step 2 of ROADMAP_TWO_WAY.md, prototype in Python:
the composition acts back on the spreading.

Klempt, Soleimani, Junker (PAMM 2023) make the front factor depend on the type
of microorganism. Here the front growth rate of the field becomes the
composition-weighted mean of species rates,

    r(x) = sum_i chi_i(x) r_i,   chi_i = phi_i / (phi_1 + phi_2),

with r_1 = (1 + d) r and r_2 = (1 - d) r (d = 0: the one-way scheme). Points
without biofilm keep the seed composition (chi = 1/2, so r unchanged there).
Everything else as in composition_spreading_fig.py: Klempt 2024 test case 4.1
field (Table 2, edge-line source, time scale 1.5, growth on both faces,
consumption g phi c), Klempt 2026 cases 3 and 6 at every node with the coupled
scheme of the ANSYS runs (amount from the field, phi_cap = 0.9, s = 0.15),
coupling step 0.1. The rates r_i are not from a paper; d is a sensitivity
parameter.

With d = 0 the field is the one-way field (r unchanged, up to round-off).

    python ansys_usermat/two_way_step2.py -> assets/fig_two_way_step2.png

d = 1/3 (main value): the cases 3 and 6 of Klempt et al. 2026 have eta_1 = 1,
eta_2 = 2, so species 1 relaxes twice as fast; r_1/r_2 = 2 gives d = 1/3. That
the front rate scales like 1/eta_i is an assumption. At coupling step 0.025:
case 3 0.704 -> 0.735, case 6 0.704 -> 0.648 (share 0.041 -> 0.046); largest
pointwise change of phi 0.32 / 0.47. Thesis figure:
    python ansys_usermat/two_way_step2.py --dt 0.025 --d 0.3333333333 --save s.npz
    python ansys_usermat/two_way_step2.py --thesis s.npz

Sensitivity, T* = 1, d = 0.5 (species 1 grows 3x as fast as species 2),
mean phi one-way -> two-way and biomass-weighted share phi_1/(phi_1+phi_2):
  coupling step   case 3                   case 6
  0.1             0.704 -> 0.735 (0.656)   0.704 -> 0.599 (0.026 -> 0.033)
  0.05            0.704 -> 0.749 (0.663)   0.704 -> 0.611 (0.034 -> 0.042)
  0.025           0.704 -> 0.750 (0.666)   0.704 -> 0.618 (0.041 -> 0.049)
The figure is made with 0.025 (python ansys_usermat/two_way_step2.py --dt 0.025);
the mean phi changes by less than 0.01 from 0.05 to 0.025. Largest pointwise
change of phi at 0.025: 0.41 (case 3), 0.55 (case 6).
Which species wins now decides how fast the biofilm spreads: in case 3 the
faster species is the majority and the front moves faster; in case 6 the
slower species takes over and the front falls behind. The composition itself
hardly changes, because it is set by the point model's own dynamics.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import composition_transport_check as T  # noqa: E402
import figstyle  # noqa: E402

OUT = HERE.parent / "assets" / "fig_two_way_step2.png"
DS = [0.0, 0.5]
CASES = {"2sp_case3": "case 3 (coexistence)", "2sp_case6": "case 6 (one species wins)"}


def run(case, d):
    B, K = T.B, T.K
    s = T.S_MACRO
    T.ms.set_case(case)
    adv = T.stepper(T.ms.ECOLOGY_CASE["theta"], T.ms.ECOLOGY_CASE["hp"])
    T.ms.set_case(None)
    phi, _, _ = K.setup("fig4_edge")
    n = phi.size
    G = np.zeros((n, 12)); started = np.zeros(n, bool)
    chi = np.full(phi.shape, 0.5)
    rec = {"t": [0.0], "phi": [phi.mean()], "chi": [0.5]}
    for k in range(1, int(round(1.0 / T.DT)) + 1):
        r_field = 100 * s * ((1 + d) * chi + (1 - d) * (1 - chi))
        K.R, K.BETA, K.K_A = r_field, 2 * s, 1e-3 * s
        _, phi, _ = B.run_setup(phi, B.nutrient_mask("edge"), 1e8, "ic", "first_order", "abs",
                                t_end=T.DT)
        p = np.minimum(phi.ravel(), T.PHI_CAP)
        new = (~started) & (p >= T.PHI_MIN)
        G[new] = T.seed_state(new.sum()); started |= new
        act = started & (p >= T.PHI_MIN)
        G[act] = T.rescale(G[act], p[act]); G[act] = T.rescale(adv(G[act]), p[act])
        c1 = np.where(started, G[:, 0] / np.maximum(G[:, 0] + G[:, 1], 1e-300), 0.5)
        chi = c1.reshape(phi.shape)
        w = np.where(started, phi.ravel(), 0.0)
        rec["t"].append(k * T.DT); rec["phi"].append(phi.mean())
        rec["chi"].append(float((w * c1).sum() / max(w.sum(), 1e-300)))
    K.R = 100.0                     # restore a scalar for other callers
    return rec, phi, chi


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--dt", type=float, default=T.DT, help="coupling step (default 0.1)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--save", help="write the results to this .npz")
    ap.add_argument("--n", type=int, default=21, help="grid nodes per edge (21 = the paper's)")
    ap.add_argument("--d", type=float, default=0.5, help="rate contrast d (1/3: from eta_1 = 1, eta_2 = 2)")
    ap.add_argument("--thesis", help="only draw the thesis figure from this .npz")
    args = ap.parse_args(argv)
    if args.thesis:
        return thesis_figure(args.thesis)
    T.DT = args.dt
    T.K.set_grid(args.n)
    DS[1] = args.d
    figstyle.apply(size=11)
    res = {(c, d): run(c, d) for c in CASES for d in DS}
    if args.save:
        np.savez(args.save, **{f"{c}|{d}|{k}": np.asarray(v) for (c, d), (rec, phi, chi) in res.items()
                               for k, v in (("t", rec["t"]), ("phi_mean", rec["phi"]), ("chi_mean", rec["chi"]),
                                            ("phi", phi), ("chi", chi))})
    for c in CASES:
        for d in DS:
            rec = res[(c, d)][0]
            print(f"{c} d={d}: mean phi at T*=1 {rec['phi'][-1]:.4f}, biomass-weighted share {rec['chi'][-1]:.3f}")
        dphi = np.abs(res[(c, DS[1])][1] - res[(c, DS[0])][1]).max()
        print(f"   largest change of phi at T*=1 for d={DS[1]}: {dphi:.3f}")
    cmap, norm = figstyle.klempt_cmap()
    K = T.K
    kz = K.N // 2
    fig, ax = plt.subplots(2, 3, figsize=(13, 8.4))
    for r, c in enumerate(CASES):
        for col, d in enumerate(DS):
            _, phi, chi = res[(c, d)]
            a = ax[r, col]
            a.imshow(np.where(phi >= 0.05, chi, np.nan)[:, :, kz].T, origin="lower",
                     extent=(0, K.L, 0, K.L), cmap=cmap, norm=norm)
            a.contour(np.linspace(0, K.L, K.N), np.linspace(0, K.L, K.N), phi[:, :, kz].T, [0.5],
                      colors="black", linewidths=1)
            a.set_facecolor("#d9dde2"); a.grid(False)
            a.set_title(("one-way ($d=0$)" if d == 0 else f"two-way ($d={d:g}$)"), fontsize=11)
            a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
        ax[r, 0].set_ylabel(CASES[c] + "\n$y$ [$\\mu$m]")
        b = ax[r, 2]
        for d, ls in zip(DS, ("-", "--")):
            rec = res[(c, d)][0]
            b.plot(rec["t"], rec["phi"], ls, color="C0", label=f"mean $\\phi$, $d={d:g}$")
            b.plot(rec["t"], rec["chi"], ls, color="C1", label=f"share $\\phi_1/(\\phi_1+\\phi_2)$, $d={d:g}$")
        b.set_xlabel("$T^*$"); b.legend(fontsize=8, frameon=False)
    for a in ax[1, :2]:
        a.set_xlabel(r"$x$ [$\mu$m]")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax[:, :2], label=r"$\phi_1/(\phi_1+\phi_2)$", shrink=0.8, pad=0.02)
    fig.suptitle("Two-way coupling, step 2 (Python): front growth rate $r=\\sum_i\\chi_i r_i$, "
                 "$r_{1,2}=(1\\pm d)\\,r$ (not from a paper);\nKlempt 2024 test case 4.1 field, "
                 "mid-plane $z=10\\ \\mu$m, black: $\\phi=0.5$, grey: no biofilm", fontsize=11.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


THESIS_OUT = HERE.parent / "assets" / "fig_two_way_step2_thesis.png"


def thesis_figure(npz):
    """Thesis version: phi at T* = 1 on the mid-plane for case 6 without and with
    the feedback, and the mean phi of both cases over time. No title (the
    caption carries it)."""
    figstyle.apply(size=11)
    z = np.load(npz)
    dk = next(k.split("|")[1] for k in z.files if k.startswith("2sp_case6|") and k.split("|")[1] != "0.0")
    dlab = "1/3" if abs(float(dk) - 1 / 3) < 1e-6 else f"{float(dk):g}"
    K = T.K
    kz = K.N // 2
    cmap, norm = figstyle.klempt_cmap()
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.9), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.15, top=0.9, wspace=0.35)
    for a, d, ttl in ((ax[0], "0.0", "one-way"), (ax[1], dk, f"with feedback, $d={dlab}$")):
        phi = z[f"2sp_case6|{d}|phi"]
        im = a.imshow(phi[:, :, kz].T, origin="lower", extent=(0, K.L, 0, K.L), cmap=cmap, norm=norm,
                      interpolation="bilinear")
        a.set_title(f"case 6, {ttl}", fontsize=11)
        a.set_xlabel(r"$x$ [$\mu$m]"); a.grid(False)
        a.set_xticks([0, 10, 20]); a.set_yticks([0, 10, 20])
    ax[0].set_ylabel(r"$y$ [$\mu$m]"); ax[1].set_yticklabels([])
    cb = fig.colorbar(im, ax=ax[:2], shrink=0.9, pad=0.02)
    cb.set_label(r"$\phi$ at $T^*=1$")
    b = ax[2]
    for c, col, lab in (("2sp_case3", "C0", "case 3"), ("2sp_case6", "C3", "case 6")):
        b.plot(z[f"{c}|0.0|t"], z[f"{c}|0.0|phi_mean"], color="0.4", lw=1.2,
               label="one-way" if c == "2sp_case3" else None)
        b.plot(z[f"{c}|{dk}|t"], z[f"{c}|{dk}|phi_mean"], "--", color=col, lw=1.8, label=f"{lab}, with feedback")
    b.set_xlabel(r"$T^*$"); b.set_ylabel(r"mean $\phi$")
    b.legend(frameon=False, fontsize=9)
    fig.savefig(THESIS_OUT, dpi=300)
    print("wrote", THESIS_OUT)


if __name__ == "__main__":
    main()
