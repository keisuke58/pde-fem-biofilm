#!/usr/bin/env python3
"""Why the four clinical conditions give almost the same stress.

The condition enters the mechanics only through alpha, and alpha only through
the living total phi_tot = sum_i phi_i psi_i (alpha_dot = k_alpha phi_tot,
ecology_jax.living_fraction_total). Which species make up that total never
reaches the stress. The CLSM compositions are fractions -- sum_i phi_i = 1 for
every condition -- so all four start from the same total and relax to the same
value, however different their make-up.

The totals part only during the early transient, mostly through the
vitalities psi_i: the alpha spread (converged; 200 and 800 steps per segment
agree) is 0.56 % at t = 1e-4, peaks at 2.4 % near t = 1e-2, and falls to
0.12 % at t = 1 and 0.003 % at t = 50. Printed by this script.

Left: the measured initial compositions (very different). Right: phi_tot(t)
for each condition, from the same seeds and ecology parameters as the
4-region all-real cylinder deck.

Second finding, printed below the plot: alpha at the deck's TIME = 1e-4 as a
function of the number of ecology steps. ANSYS takes 10 (NSUBST,10, one
ecology step per increment) and reproduces the n=10 column exactly; the
converged values differ from it by ~1 %, more than the spread between the
conditions, and in a different order. The spread itself (~0.5 %) survives
refinement; the ranking of the conditions does not.

    python ansys_usermat/apdl/plot_condition_total_phi.py -o assets/condition_total_phi.png
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "coupling"))
import ecology_deck_reference as R  # noqa: E402
import ecology_jax as E  # noqa: E402

DECK = HERE / "t_growth_cylinder_ecology_4region_all_real_clsm.dat"
CONDS = {2: "CS", 3: "CH", 4: "DS", 5: "DH"}
SPECIES = ["S. oralis", "A. naeslundii", "Veillonella", "F. nucleatum",
           "P. gingivalis"]
# Categorical slots 1-5 in fixed order (dataviz reference palette, light).
SPECIES_COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
DECK_TIME = 1e-4


def seeds():
    lines = DECK.read_text(errors="replace").splitlines()
    prop = R.table(lines, "USER", 2)
    k_alpha, theta = prop[8], prop[9:29]
    g0 = {}
    for mat, name in CONDS.items():
        s = np.zeros(26)
        v = R.table(lines, "STATE", mat)
        s[:len(v)] = v
        g0[name] = s[14:26]
    return k_alpha, theta, g0


def trajectory(g0, theta, t_edges, per_seg=200):
    """(phi_tot, int_0^t phi_tot) at each edge, integrating each log-spaced
    segment in per_seg equal steps."""
    g, t_prev, out, cum = np.asarray(g0), 0.0, [], []
    acc = 0.0
    for t in t_edges:
        g, pi = E.ecology_substeps(g, theta, t - t_prev, per_seg)
        acc += pi
        out.append(E.living_fraction_total(g))
        cum.append(acc)
        t_prev = t
    return np.array(out), np.array(cum)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/condition_total_phi.png"))
    ap.add_argument("--per-seg", type=int, default=200,
                    help="ecology steps per log-spaced time segment "
                         "(convergence check: compare 200 with 800)")
    a = ap.parse_args(argv)

    k_alpha, theta, g0 = seeds()
    t = np.logspace(-7, np.log10(50.0), 70)
    res = {c: trajectory(g, theta, t, a.per_seg) for c, g in g0.items()}
    traj = {c: r[0] for c, r in res.items()}

    print(f"alpha spread between the conditions vs growth time "
          f"({a.per_seg} steps per log segment):")
    for T in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 50.0):
        j = int(np.argmin(abs(np.log(t / T))))
        al = {c: k_alpha * r[1][j] for c, r in res.items()}
        lo, hi = min(al, key=al.get), max(al, key=al.get)
        print(f"  t={t[j]:9.3e}: spread {(al[hi] / al[lo] - 1) * 100:6.3f} %"
              f"  ({lo} lowest, {hi} highest)")

    print(f"alpha at TIME={DECK_TIME:g} vs number of ecology steps "
          f"(ANSYS uses 10):")
    ns = [1, 10, 100, 1000, 10000]
    print("cond " + "".join(f"{'n=' + str(n):>13}" for n in ns))
    conv = {}
    for c, g in g0.items():
        row = [k_alpha * E.ecology_substeps(g, theta, DECK_TIME, n)[1] for n in ns]
        conv[c] = row
        print(f"{c:<5}" + "".join(f"{v:13.5e}" for v in row))
    for j, n in ((1, 10), (4, 10000)):
        v = [conv[c][j] for c in conv]
        order = sorted(conv, key=lambda c: conv[c][j])
        print(f"n={n}: spread {(max(v) / min(v) - 1) * 100:.2f} %, "
              f"order low->high {' < '.join(order)}")
    print(f"phi_tot at t={t[-1]:g}: " + ", ".join(
        f"{c} {traj[c][-1]:.4f}" for c in traj))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"axes.edgecolor": "0.6", "axes.labelcolor": "0.2",
                         "xtick.color": "0.35", "ytick.color": "0.35"})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.0),
                                   gridspec_kw={"width_ratios": [1, 1.5]})

    names = list(g0)
    for i, c in enumerate(names):
        bottom = 0.0
        for j in range(5):
            h = g0[c][j]
            if h > 0:
                ax1.bar(i, h, bottom=bottom, width=0.62, color=SPECIES_COL[j],
                        edgecolor="white", linewidth=1.5,
                        label=SPECIES[j] if i == 0 else None)
            bottom += h
    ax1.set_xticks(range(len(names)))
    ax1.set_xticklabels(names)
    ax1.set_ylim(0, 1.0)
    ax1.set_ylabel(r"initial fraction $\varphi_i$")
    ax1.set_title("Composition: very different", fontsize=11)
    ax1.legend(fontsize=7, frameon=False, loc="upper left",
               bbox_to_anchor=(1.0, 1.0), handlelength=1.0)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    styles = ["-", "--", ":", "-."]
    for (c, y), ls in zip(traj.items(), styles):
        ax2.plot(t, y, ls, lw=2, color="0.25", label=c)
    ax2.axvline(DECK_TIME, color="0.6", lw=1)
    ax2.annotate("cylinder decks\n(TIME = 1e-4)", (DECK_TIME, 0.30),
                 xytext=(6, 0), textcoords="offset points", fontsize=8,
                 color="0.35")
    ax2.set_xscale("log")
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel("time $t$")
    ax2.set_ylabel(r"$\varphi_{tot}=\sum_i \varphi_i\psi_i$  (drives $\dot\alpha$)")
    ax2.set_title(r"What $\alpha$ sees: the total -- the same", fontsize=11)
    ax2.legend(fontsize=8, frameon=False, loc="lower left")
    ax2.grid(alpha=0.25)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=200)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
