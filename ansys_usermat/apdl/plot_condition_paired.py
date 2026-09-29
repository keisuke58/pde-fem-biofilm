#!/usr/bin/env python3
"""Separating condition from position, using two runs that already exist.

One four-region run cannot say whether the stress differences between the
clinical conditions are caused by the conditions or by where on the shell each
one happens to sit: each condition occupies one quadrant, so the two vary
together. Chapter 5 said so -- and then went one step further, reporting that
position accounted for most of the variation. That step was not supported by
the single run either, and it was wrong.

Two runs settle it, and both are already in this repository:

    t_growth_cylinder_ecology_4region_clsm.dat            A=default B=CS C=CH D=DS
    t_growth_cylinder_ecology_4region_all_real_clsm.dat   A=CS  B=CH  C=DS  D=DH

The decks are identical apart from the seeds -- same geometry, same mesh, same
material constants, same time stepping -- and the assignment is shifted by one
material slot, so CS, CH and DS each appear in two different quadrants. This
script checks that the element sets really are identical before comparing
anything, then reads the same condition at two placements.

The answer is the opposite of what the chapter claimed. Holding the condition
and Z fixed and moving only in theta changes the mean von Mises by about a
tenth of a percent, while the conditions differ from each other by about nine
percent.

One thing this does not explain, and the chapter should not pretend it does:
the stress is not monotone in the growth. CS has a larger alpha than DS and a
lower stress. The material constants are identical across the regions, so
something in the composition reaches the stress by a route other than alpha;
that is left open rather than guessed at.

    python ansys_usermat/apdl/plot_condition_paired.py -o assets/v222_conditions_paired.png
"""
import argparse
import re
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_condition_comparison import parse  # noqa: E402

# Quadrant of each material, from the decks' shared EMODIF block.
QUAD = {2: ("low", "low"), 3: ("high", "low"),
        4: ("low", "high"), 5: ("high", "high")}

# Which condition sat in which material, per run. Verified below against the
# alpha the solver actually reported, so a mislabelling here cannot pass.
RUNS = {
    "growth_result_cylinder_ecology_4region_clsm.txt":
        {2: "default", 3: "CS", 4: "CH", 5: "DS"},
    "growth_result_cylinder_ecology_4region_all_real_clsm.txt":
        {2: "CS", 3: "CH", 4: "DS", 5: "DH"},
}
ALPHA = {"CS": 4.6145e-3, "CH": 4.6920e-3, "DS": 4.3557e-3,
         "DH": 4.8186e-3, "default": 4.5687e-3}
ORDER = ["CS", "DS", "default", "CH", "DH"]


def read_alpha(path):
    """alpha per material, from the SVAR(10) listing. Uniform within a region."""
    txt = Path(path).read_text(errors="replace")
    out, parts = {}, re.split(r"IN RANGE\s+(\d+) TO\s+\d+", txt)
    for i in range(1, len(parts), 2):
        seg = parts[i + 1]
        m = re.search(r"NODE\s+10\s*\n", seg)
        if not m:
            continue
        vals = {float(v) for _, v in
                re.findall(r"\s+(\d+)\s+([-+0-9.E]+)\s*\n", seg[m.end():m.end() + 4000])}
        if len(vals) == 1:
            out[int(parts[i])] = vals.pop()
    return out


def collect(here):
    """[(condition, theta, Z, n, mean SEQV)], both runs, checked as it goes."""
    sets, rows = None, []
    for fn, assign in RUNS.items():
        seqv, mats = parse(here / fn)
        if sets is None:
            sets = mats
        elif any(sets[m] != mats[m] for m in sets):
            raise SystemExit("the two runs do not cover the same elements per "
                             "material; they are not comparable")
        alpha = read_alpha(here / fn)
        for m in sorted(mats):
            cond = assign[m]
            if m in alpha and abs(alpha[m] - ALPHA[cond]) > 1e-9:
                raise SystemExit(
                    f"{fn} material {m}: the deck's own alpha {alpha[m]:.6e} is "
                    f"not {cond}'s {ALPHA[cond]:.6e} -- the assignment table "
                    f"in this script is wrong, fix it before trusting anything")
            rows.append((cond, QUAD[m][0], QUAD[m][1], len(mats[m]),
                         st.mean(seqv[e] for e in mats[m])))
    return rows


def theta_contrasts(rows):
    """Same condition, same Z, theta differing: the only clean position test."""
    out = []
    for c in {r[0] for r in rows}:
        rs = [r for r in rows if r[0] == c]
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                if rs[i][2] == rs[j][2] and rs[i][1] != rs[j][1]:
                    lo = rs[i] if rs[i][1] == "low" else rs[j]
                    hi = rs[j] if rs[i][1] == "low" else rs[i]
                    out.append((c, lo[2], hi[4] / lo[4]))
    return sorted(out)


def plot(rows, out, caption=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conds = [c for c in ORDER if any(r[0] == c for r in rows)]
    # Without the caption the figure can be shorter: on a slide the block is
    # unreadable at the size the frame allows, and the speaker says it instead.
    fig, ax = plt.subplots(figsize=(9.6, 4.6) if caption else (9.6, 3.6))

    for i, c in enumerate(conds):
        rs = sorted((r for r in rows if r[0] == c), key=lambda r: r[4])
        ys = [r[4] * 1e6 for r in rs]
        col = "0.45" if c == "default" else "#c0392b"
        if len(ys) > 1:
            ax.plot([i, i], [min(ys), max(ys)], lw=2.2, color=col, zorder=2,
                    solid_capstyle="round")
        # Where the two placements nearly coincide -- which is the finding --
        # their labels would sit on top of each other, so stagger them.
        tight = len(ys) > 1 and (max(ys) - min(ys)) < 0.05 * max(ys)
        for k, (r, y) in enumerate(zip(rs, ys)):
            ax.plot(i, y, marker="o", ms=8, color=col, zorder=3)
            dy = (9 if k else -9) if tight else 0
            ax.annotate(f"$\\theta$ {r[1]}, $Z$ {r[2]}", (i, y),
                        xytext=(11, dy), textcoords="offset points",
                        fontsize=7.4, va="center", color="0.35")

    vals = [r[4] for r in rows if r[0] != "default"]
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds, fontsize=10)
    ax.set_xlim(-0.5, len(conds) - 0.15)
    allv = [r[4] * 1e6 for r in rows]
    pad = 0.09 * (max(allv) - min(allv))     # room for the staggered labels
    ax.set_ylim(min(allv) - pad, max(allv) + pad)
    ax.set_ylabel(r"mean von Mises  [$\times 10^{-6}$]")
    ax.set_title("The same condition, put in two different places, gives the "
                 "same stress", fontsize=11)
    ax.grid(alpha=0.25, axis="y")

    if caption:
        tc = theta_contrasts(rows)
        txt = "; ".join(f"{c} at $Z$ {z}: {(r - 1) * 100:+.2f}%" for c, z, r in tc)
        fig.text(0.5, 0.02,
                 f"Each vertical bar joins one condition's two placements across "
                 f"the two runs (identical mesh, materials and stepping).\n"
                 f"Moving in $\\theta$ at fixed condition and $Z$ -- {txt} -- "
                 f"against a spread of {(max(vals) / min(vals) - 1) * 100:.1f}% "
                 f"between conditions.\n"
                 f"Position is worth about a tenth of a percent here; the "
                 f"conditions are what separate.",
                 ha="center", va="bottom", fontsize=8, color="0.35")
    fig.tight_layout(rect=(0, 0.16 if caption else 0.0, 1, 1))
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_conditions_paired.png"))
    ap.add_argument("--no-caption", action="store_true",
                    help="omit the explanatory block; for slides, where it is "
                         "unreadable at the size the frame allows")
    a = ap.parse_args(argv)

    rows = collect(here)
    print(f"{'cond':<9}{'theta':<7}{'Z':<6}{'n':>3}  mean SEQV")
    for r in sorted(rows):
        print(f"{r[0]:<9}{r[1]:<7}{r[2]:<6}{r[3]:>3}  {r[4]:.5e}")
    print("\nclean theta contrasts (same condition, same Z):")
    for c, z, ratio in theta_contrasts(rows):
        print(f"  {c} at Z {z}: high/low = {ratio:.5f}  ({(ratio - 1) * 100:+.2f}%)")
    vals = [r[4] for r in rows if r[0] != "default"]
    print(f"\nspread between conditions: {max(vals) / min(vals):.4f} "
          f"({(max(vals) / min(vals) - 1) * 100:.1f}%)")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(rows, a.out, caption=not a.no_caption)
    return 0


if __name__ == "__main__":
    sys.exit(main())
