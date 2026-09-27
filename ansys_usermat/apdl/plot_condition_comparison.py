#!/usr/bin/env python3
"""The four clinical conditions, as ANSYS solved them, and what they cannot show.

`t_growth_cylinder_ecology_4region_all_real_clsm.dat` seeds four regions of one
curved shell with the four real Day-1 CLSM compositions -- CS, CH, DS, DH --
drives each through the 0D Hamilton ecology at the Gauss point, and reports the
resulting growth and von Mises stress. It is the run Chapter 5's results section
was waiting for, and this reads its listing.

It also reads the limit the run has, which matters more than the numbers.
The four regions are the four quadrants of a 2x2 (theta x Z) split of one
shell, one condition per quadrant, so **condition is perfectly confounded with
position**: every difference between two conditions is also a difference
between two places on the geometry, and no amount of arithmetic separates them.

That is not a hypothetical worry here. The measured stress ordering does not
follow the growth ordering -- DS has the lowest alpha of the four and the
second-highest mean stress -- and splitting the same numbers by position
instead of by condition explains them better: the high-theta half is 7.3 %
above the low-theta half, against a total spread of 8.9 % across all four
cells. Position is doing most of the work.

So the figure reports the growth per condition, which is sound, and shows the
stress in the physical 2x2 layout rather than as a condition ranking, which
would invite a comparison the design cannot support. Getting a real condition
comparison needs the same region solved four times with different seeds, not
four regions solved once -- one ANSYS session, once the machine is available.

    python ansys_usermat/apdl/plot_condition_comparison.py -o assets/v222_conditions.png
"""
import argparse
import re
import statistics as st
import sys
from pathlib import Path

# mat -> (condition, theta half, Z half), read off the deck's EMODIF block
LAYOUT = {2: ("CS", "Commensal/Static", "Comm./Static", 0, 0),
          3: ("CH", "Commensal/HOBIC", "Comm./HOBIC", 1, 0),
          4: ("DS", "Dysbiotic/Static", "Dysb./Static", 0, 1),
          5: ("DH", "Dysbiotic/HOBIC", "Dysb./HOBIC", 1, 1)}

# alpha after 10 chained substeps, from ecology_4region_reference_all_real_clsm.py
# -- the independent reference the ANSYS run was verified against, not a
# transcription of the run's own output.
ALPHA = {2: 4.6144769029e-03, 3: 4.6920293146e-03,
         4: 4.3556778678e-03, 5: 4.8186255156e-03}


def parse(path):
    """Element von Mises, and the element set of each material.

    SEQV is taken only from inside ELEMENT TABLE LISTING sections: the SVAR
    listings further down are also `number number` rows and would otherwise be
    read as stresses. The parse is checked against the listing's own MINIMUM
    and MAXIMUM lines before anything is plotted.
    """
    txt = Path(path).read_text(errors="replace")
    seqv = {}
    for sec in re.findall(r"POST1 ELEMENT TABLE LISTING.*?"
                          r"(?=POST1|MINIMUM VALUES|\Z)", txt, re.S):
        for m in re.finditer(r"^\s+(\d+)\s+(-?\d*\.\d+E[+-]\d+)\s*$", sec, re.M):
            seqv[int(m.group(1))] = float(m.group(2))
    if not seqv:
        raise SystemExit(f"no element table found in {path}")

    lo = re.search(r"MINIMUM VALUES\s*\n\s*ELEM\s+(\d+)\s*\n\s*VALUE\s+"
                   r"(-?\d*\.\d+E[+-]\d+)", txt)
    hi = re.search(r"MAXIMUM VALUES\s*\n\s*ELEM\s+(\d+)\s*\n\s*VALUE\s+"
                   r"(-?\d*\.\d+E[+-]\d+)", txt)
    for tag, m in (("min", lo), ("max", hi)):
        if m and abs(seqv.get(int(m.group(1)), 0) - float(m.group(2))) > 0:
            raise SystemExit(f"parse disagrees with the listing's own {tag} "
                             f"at element {m.group(1)}")

    mats, parts = {}, re.split(r"IN RANGE\s+(\d+) TO\s+\d+", txt)
    for i in range(1, len(parts), 2):
        els = sorted(set(int(x) for x in
                         re.findall(r"ELEMENT=\s*(\d+)\s+SOLID185", parts[i + 1])))
        if els:
            mats[int(parts[i])] = els
    return seqv, mats


def plot(seqv, mats, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mean = {m: st.mean(seqv[e] for e in mats[m]) for m in mats}
    sub = sorted(set(seqv) - set().union(*mats.values()))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.5),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    # ---- left: growth per condition, the part that is a clean comparison ----
    keys = [m for m in (2, 3, 4, 5) if m in mats]   # CS, CH, DS, DH
    vals = [ALPHA[m] for m in keys]
    axL.bar(range(len(keys)), vals, width=0.6, color="#5b7fa6",
            edgecolor="white", lw=1.5)
    for i, (m, v) in enumerate(zip(keys, vals)):
        axL.text(i, v * 1.005, f"{v:.4e}", ha="center", va="bottom", fontsize=8,
                 color="0.25")
    axL.set_xticks(range(len(keys)))
    axL.set_xticklabels([f"{LAYOUT[m][0]}\n{LAYOUT[m][2]}" for m in keys],
                        fontsize=8.5)
    axL.set_ylabel(r"growth $\alpha$ after 10 substeps")
    axL.set_ylim(min(vals) * 0.93, max(vals) * 1.04)
    axL.set_title("Four real Day-1 CLSM compositions drive four\n"
                  "different growth values", fontsize=10.5)
    axL.grid(alpha=0.25, axis="y")
    # Anchored to the figure, not the axes: axes-fraction text moves with the
    # axes when tight_layout shrinks them, and lands back on the tick labels.
    fig.text(0.27, 0.055,
             f"spread {max(vals) / min(vals):.3f}x across the four. Each value "
             f"matches an independent\nreference implementation to full printed "
             f"precision.",
             ha="center", va="bottom", fontsize=7.8, color="0.35")

    # ---- right: stress in the physical layout, not as a condition ranking ----
    lo, hi = min(mean.values()), max(mean.values())
    for m in mats:
        _, _, _, th, z = LAYOUT[m]
        f = (mean[m] - lo) / (hi - lo) if hi > lo else 0.5
        axR.add_patch(plt.Rectangle((th, z), 1, 1,
                                    fc=plt.get_cmap("YlOrRd")(0.25 + 0.55 * f),
                                    ec="white", lw=2.5))
        axR.text(th + 0.5, z + 0.62, LAYOUT[m][0], ha="center", va="center",
                 fontsize=13, weight="bold", color="0.15")
        axR.text(th + 0.5, z + 0.38, f"{mean[m]:.3e}", ha="center", va="center",
                 fontsize=8.5, color="0.2")
    axR.set_xlim(-0.05, 2.05)
    axR.set_ylim(-0.05, 2.05)
    axR.set_xticks([0.5, 1.5])
    axR.set_xticklabels(["low " + r"$\theta$", "high " + r"$\theta$"])
    axR.set_yticks([0.5, 1.5])
    axR.set_yticklabels(["low $Z$", "high $Z$"])
    axR.tick_params(length=0)
    for s in axR.spines.values():
        s.set_visible(False)
    axR.set_title("...but the same stresses split better by position\n"
                  "than by condition", fontsize=10.5)

    def half(idx, v):
        return st.mean(mean[m] for m in mats if LAYOUT[m][3 + idx] == v)
    th_r = half(0, 1) / half(0, 0)
    z_r = half(1, 1) / half(1, 0)
    fig.text(0.77, 0.025,
             f"mean von Mises per quadrant. high-$\\theta$ / low-$\\theta$ = "
             f"{th_r:.3f},  high-$Z$ / low-$Z$ = {z_r:.3f},\n"
             f"against a total spread of {hi / lo:.3f} across all four. One "
             f"condition per quadrant means\ncondition and position cannot be "
             f"separated -- this is not a condition comparison yet.",
             ha="center", va="bottom", fontsize=7.8, color="0.35")

    fig.suptitle("Four clinical conditions on one curved shell, ANSYS 2022 R2 "
                 f"({len(seqv)} elements)", fontsize=11)
    fig.tight_layout(rect=(0, 0.20, 1, 0.93))
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    return mean, sub


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--listing", type=Path, default=here /
                    "growth_result_cylinder_ecology_4region_all_real_clsm.txt")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_conditions.png"))
    a = ap.parse_args(argv)

    seqv, mats = parse(a.listing)
    missing = set(LAYOUT) - set(mats)
    if missing:
        raise SystemExit(f"listing has no elements for material(s) {missing}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    mean, sub = plot(seqv, mats, a.out)

    print(f"{len(seqv)} elements; {len(sub)} of them substrate (material 1)")
    for m in sorted(mats):
        print(f"  {LAYOUT[m][0]}  {LAYOUT[m][1]:<20} n={len(mats[m]):>2}  "
              f"alpha={ALPHA[m]:.4e}  mean SEQV={mean[m]:.4e}")
    if sub:
        sm = st.mean(seqv[e] for e in sub)
        print(f"  substrate            n={len(sub):>2}  alpha=0            "
              f"mean SEQV={sm:.4e}  "
              f"({st.mean(mean.values()) / sm:.1f}x below the growth layer)")
    order_a = [LAYOUT[m][0] for m in sorted(mats, key=lambda m: ALPHA[m])]
    order_s = [LAYOUT[m][0] for m in sorted(mats, key=lambda m: mean[m])]
    print(f"\nalpha order  {' < '.join(order_a)}")
    print(f"stress order {' < '.join(order_s)}")
    if order_a != order_s:
        print("The two orders differ, which is the confound: condition and "
              "position vary together in this deck.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
