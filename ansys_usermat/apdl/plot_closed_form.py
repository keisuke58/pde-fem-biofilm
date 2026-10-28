#!/usr/bin/env python3
"""The two growth verification cases, as ANSYS 2022 R2 actually solved them.

`t_growth_free.dat` and `t_growth_constrained.dat` are the same single element,
the same material (C10=2e-4, D1=5e3, eta=0, alpha=0.05) and the same growth
drive. Only the boundary conditions differ, and they differ at the two
extremes:

  free         only the six rigid-body modes are removed, so the element grows
               into the shape Fg wants. Fe = I, and the stress is EXACTLY zero
               for any material and any alpha.
  constrained  every node is fixed, F = I, so Fe = I/(1+alpha). The elastic
               state is isotropic, its deviator vanishes, and only the
               volumetric term survives.

The pair is not redundant. Full constraint can mask a sign error in
Fe = F Fg^-1 -- the wrong sign cancels when the element cannot move -- and the
free case is what catches it. This is why both are run.

The right-hand panel is the more interesting one. The constrained answer the
solver returns is NOT the ideal continuum value: it is the ideal value plus the
spurious spherical term of DEVIATOR_SCALING_FINDING.md, which is predicted in
closed form rather than merely observed. Their sum reproduces the solver to
every digit MAPDL prints. That term is 47% of the reported stress here, which
is worth showing plainly -- it is a pressure error, so the von Mises stress
this work reports is unaffected, but the hydrostatic number is not the ideal
model's number and the figure should not pretend otherwise.

Nothing here imports the implementation; the predictions come from
closed_form_reference.py, which is derived from the continuum statement.

    python ansys_usermat/apdl/plot_closed_form.py -o assets/v222_closed_form.png
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from closed_form_reference import constrained_stress, spurious_term  # noqa: E402

ALPHA, D1, C10 = 0.05, 5.0e3, 0.2e-3
COMPONENTS = ["SX", "SY", "SZ", "SXY", "SYZ", "SXZ"]

# A stress row is a node number followed by six numbers. MAPDL prints them in
# fixed-width columns with no separating space when a value is negative
# ("-0.19233E-010-0.19233E-010"), so the row cannot be split on whitespace --
# it has to be scanned for number tokens.
_NUM = re.compile(r"[-+]?\d*\.\d+E[-+]\d+|[-+]?\d+\.\d+|[-+]?\d+")
_HDR = re.compile(r"NODE\s+SX\s+SY\s+SZ")


def parse_stress(path):
    """Returns {node: [SX, SY, SZ, SXY, SYZ, SXZ]} from the first stress table.

    Stops at the first blank line after the table so the SVAR listings further
    down the file, which have a different column layout, cannot be misread as
    stress rows.
    """
    rows, active = {}, False
    for line in Path(path).read_text(errors="replace").splitlines():
        if _HDR.search(line):
            active = True
            continue
        if not active:
            continue
        toks = _NUM.findall(line)
        if len(toks) == 7:
            rows[int(float(toks[0]))] = [float(t) for t in toks[1:]]
        elif rows and not line.strip():
            break
    if not rows:
        raise SystemExit(f"no stress table found in {path}")
    return rows


def plot(free, cons, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import figstyle
    figstyle.apply(8.5)          # drawn at the thesis text width (6.3 in): text prints at 7.5-9.5 pt
    import numpy as np

    ideal = constrained_stress(ALPHA, D1)[0]
    spur = spurious_term(ALPHA, C10)[0]
    pred = ideal + spur

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.3, 2.75),
                                   gridspec_kw={"width_ratios": [1.2, 1]})

    # ---- left: every stress component of both solves, on one log axis ----
    for k, (name, rows, col) in enumerate(
            [("free growth", free, "#1e8449"),
             ("constrained growth", cons, "#c0392b")]):
        allv = [abs(v) for r in rows.values() for v in r]
        vals = [v for v in allv if v > 0.0]
        x = np.full(len(vals), k) + np.linspace(-0.17, 0.17, len(vals))
        axL.scatter(x, vals, s=9, color=col, alpha=0.75, zorder=3, linewidths=0)

    axL.axhline(abs(pred), ls="--", lw=1.0, color="#c0392b", zorder=2)
    axL.text(1.42, abs(pred) * 0.3, "closed form,\nconstrained", ha="right", va="top", fontsize=7.5,
             color="#c0392b")
    axL.text(0.0, 1e-9, "closed form: 0\n(solver noise)", ha="center", va="bottom", fontsize=7.5,
             color="#1e8449")
    axL.set_yscale("log")
    axL.set_ylim(1e-23, 1e-3)
    axL.set_yticks([1e-21, 1e-15, 1e-9, 1e-3])
    axL.set_xlim(-0.45, 1.45)
    axL.set_xticks([0, 1])
    axL.set_xticklabels(["free", "constrained"])
    axL.set_ylabel(r"$|\sigma_{ij}|$")
    axL.set_title("(a) all stress components", fontsize=9)
    axL.grid(alpha=0.25, axis="y", which="major")
    for lab, col in zip(axL.get_xticklabels(), ["#1e8449", "#c0392b"]):
        lab.set_color(col)

    # ---- right: what the constrained number is actually made of ----
    meas = np.mean([r[0] for r in cons.values()])
    bars = [("continuum term", ideal, "#5b7fa6"), ("spherical term", spur, "#d68910")]
    bottom = 0.0
    for label, v, col in bars:
        axR.bar(0, -v * 1e4, bottom=-bottom * 1e4, width=0.5, color=col, edgecolor="white", lw=1.0,
                label=label)
        bottom += v
    axR.bar(1, -meas * 1e4, width=0.5, color="0.35", edgecolor="white", lw=1.0, label="ANSYS")
    axR.axhline(-pred * 1e4, ls="--", lw=1.0, color="0.2", zorder=4)
    axR.set_xticks([0, 1])
    axR.set_xticklabels(["closed form", "ANSYS"])
    axR.set_ylabel(r"$-\sigma_{xx}\ [10^{-4}]$")
    axR.set_title("(b) constrained value", fontsize=9)
    axR.set_ylim(0, -pred * 1e4 * 1.6)
    axR.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axR.grid(alpha=0.25, axis="y")
    axR.legend(fontsize=7.5, frameon=False, loc="upper center", ncol=3, handlelength=1.0, columnspacing=0.8)

    fig.tight_layout(w_pad=1.5)
    fig.savefig(out, dpi=250); fig.savefig(Path(out).with_suffix(".pdf"))   # PDF: vector, used by the thesis
    print(f"closed form {pred:.6e}, ANSYS {meas:.6e}, spherical term {100 * spur / pred:.0f} %")
    print(f"wrote {out}")


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--free", type=Path, default=here / "growth_free_result.txt")
    ap.add_argument("--constrained", type=Path,
                    default=here / "growth_result.txt")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_closed_form.png"))
    a = ap.parse_args(argv)

    free, cons = parse_stress(a.free), parse_stress(a.constrained)
    ideal = constrained_stress(ALPHA, D1)[0]
    spur = spurious_term(ALPHA, C10)[0]
    pred = ideal + spur
    meas = sum(r[0] for r in cons.values()) / len(cons)
    print(f"free:        {len(free)} nodes, max |sigma| = "
          f"{max(abs(v) for r in free.values() for v in r):.3e} (closed form 0)")
    print(f"constrained: {len(cons)} nodes, sigma_xx = {meas:.6e}")
    print(f"             ideal {ideal:.6e} + spurious {spur:.6e} = {pred:.6e}")
    rel = abs(meas - pred) / abs(pred)
    print(f"             relative difference {rel:.2e}")
    if rel > 1e-4:
        print("NOTE: measured and predicted differ by more than the printed "
              "precision of the listing -- re-check before trusting the caption")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(free, cons, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
