#!/usr/bin/env python3
"""Newton convergence from the v222 solver logs.

Chapter 5 argues the finite-difference tangent is accurate enough (1.1e-6
against exact AD over 648 parameter combinations). That is a statement about
the tangent in isolation; what a reader actually wants to know is whether a
real solve converges with it, because a subtly wrong tangent does not give a
wrong answer -- it gives an iteration that will not converge.

`out_cyl_wrapper.txt` -- the 12240-element run -- answers it, and the answer
has two halves, which is why the figure has two panels:

  left   within each substep the force residual drops below the convergence
         criterion in one to four equilibrium iterations. The drop is not
         strictly monotone -- substep 1 rises once, from 1.70e-6 to 1.96e-6,
         between iterations 2 and 3 -- which is ordinary for a Newton step
         under an inexact tangent and is not smoothed away here.
  right  AUTOTS was on, so the solver was free to bisect if the tangent had
         been poor. It never did: it GREW the increment at every opportunity,
         0.5 -> 0.75 -> 1.125 -> 1.6875, finishing the load step in 6 substeps
         and 13 cumulative iterations. The last bar is short only because the
         load step ends at TIME=5; the increment ANSYS would have taken had
         there been more to solve was still growing.

Residuals are plotted per substep rather than concatenated in file order. The
concatenated view has a sawtooth -- each new substep starts from a fresh, large
residual -- which reads as non-convergence when it is the opposite.

The single-element logs deliberately do not appear: that case is fully
constrained, so it has no free degrees of freedom and its force residual is
identically 0.000 in both builds. The two agreeing there is true and vacuous,
and plotting it would suggest evidence that is not present.

    python ansys_usermat/apdl/plot_convergence.py -o assets/v222_convergence.png
"""
import argparse
import re
import sys
from pathlib import Path

_REC = re.compile(r"FORCE CONVERGENCE VALUE\s*=\s*([0-9.EDed+-]+)\s+"
                  r"CRITERION=\s*([0-9.EDed+-]+)(\s*<<< CONVERGED)?")
_DONE = re.compile(r"LOAD STEP\s+(\d+)\s+SUBSTEP\s+(\d+)\s+COMPLETED\."
                   r"\s+CUM ITER\s*=\s*(\d+)")
_TIME = re.compile(r"\*\*\*\s*TIME\s*=\s*([0-9.EDed+-]+)\s+"
                   r"TIME INC\s*=\s*([0-9.EDed+-]+)")


def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse(path):
    """Group the convergence records by substep.

    Returns a list of substeps, each {'n', 'cum_iter', 'time', 'inc', 'rec'},
    where 'rec' is [(residual, criterion, converged), ...] in solver order --
    the first entry is the residual the substep starts from, before equilibrium
    iteration 1. A trailing group with no SUBSTEP COMPLETED line (an aborted
    run) is dropped, and the caller is told.
    """
    substeps, buf, pending = [], [], None
    for line in Path(path).read_text(errors="replace").splitlines():
        m = _REC.search(line)
        if m:
            buf.append((_f(m.group(1)), _f(m.group(2)), bool(m.group(3))))
            continue
        m = _DONE.search(line)
        if m:
            pending = {"n": int(m.group(2)), "cum_iter": int(m.group(3)),
                       "rec": buf}
            buf = []
            continue
        m = _TIME.search(line)
        if m and pending is not None:
            pending["time"], pending["inc"] = _f(m.group(1)), _f(m.group(2))
            substeps.append(pending)
            pending = None
    if not substeps:
        raise SystemExit(f"no completed substeps in {path}")
    return substeps, len(buf)


def plot(substeps, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.3),
                                   gridspec_kw={"width_ratios": [1.45, 1]})

    cmap = plt.get_cmap("viridis")
    n = len(substeps)
    crit = substeps[0]["rec"][0][1]
    for k, s in enumerate(substeps):
        res = [r[0] for r in s["rec"]]
        x = list(range(len(res)))          # 0 = start of substep
        col = cmap(0.08 + 0.78 * k / max(n - 1, 1))
        axL.semilogy(x, res, marker="o", ms=4.5, lw=1.4, color=col,
                     label=f"substep {s['n']}")
        if s["rec"][-1][2]:
            axL.plot(x[-1], res[-1], marker="o", ms=11, mfc="none",
                     mec="#1e8449", mew=1.6, zorder=5)

    axL.axhline(crit, ls="--", lw=1.3, color="#c0392b")
    axL.text(0.99, crit * 1.25, "ANSYS convergence criterion",
             transform=axL.get_yaxis_transform(), ha="right", va="bottom",
             fontsize=8.5, color="#c0392b")
    axL.set_xlabel("equilibrium iteration within the substep")
    axL.set_ylabel("force residual")
    axL.set_title("Residual falls to the criterion in every substep",
                  fontsize=10.5)
    axL.set_xticks(range(max(len(s["rec"]) for s in substeps)))
    axL.grid(alpha=0.25, which="both")
    axL.legend(fontsize=8, frameon=False, loc="lower right", ncol=2)
    axL.set_ylim(top=crit * 6)

    inc = [s["inc"] for s in substeps]
    ns = [s["n"] for s in substeps]
    axR.bar(ns, inc, width=0.62, color="#5b7fa6", edgecolor="white", lw=1.2)
    for s in substeps:
        axR.text(s["n"], s["inc"] * 1.02, f"{len(s['rec']) - 1}",
                 ha="center", va="bottom", fontsize=8, color="0.3")
    axR.set_xlabel("substep")
    axR.set_ylabel("time increment applied")
    axR.set_title("AUTOTS grew the increment; it never bisected", fontsize=10.5)
    # The final increment is short because the load step ends at TIME=5, not
    # because the solver backed off -- say so, or the last bar reads as a
    # cut-back, which is the one thing this panel exists to rule out.
    axR.text(ns[-1], inc[-1] * 1.45, "short only because\nthe load step ends",
             ha="center", va="bottom", fontsize=8, color="0.35")
    axR.set_xticks(ns)
    axR.set_ylim(0, max(inc) * 1.3)
    axR.grid(alpha=0.25, axis="y")
    axR.text(0.03, 0.95,
             "numbers above the bars are equilibrium\n"
             "iterations needed for that substep",
             transform=axR.transAxes, fontsize=8, va="top", color="0.35")

    cum = substeps[-1]["cum_iter"]
    fig.suptitle(f"Newton convergence, 12240 elements, delivered routine "
                 f"-- load step finished in {len(substeps)} substeps, "
                 f"{cum} cumulative iterations", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", type=Path, default=here / "out.txt")
    ap.add_argument("--wrapper", type=Path, default=here / "out_wrapper.txt")
    ap.add_argument("--multi", type=Path, default=here / "out_cyl_wrapper.txt")
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_convergence.png"))
    a = ap.parse_args(argv)

    for p in (a.core, a.wrapper):
        txt = Path(p).read_text(errors="replace")
        vals = {_f(m.group(1)) for m in _REC.finditer(txt)}
        if vals - {0.0}:
            print(f"NOTE: {p.name} has a non-zero residual {vals} -- that case "
                  f"was assumed fully constrained; re-check the caption")
        else:
            print(f"{p.name}: residuals {vals} (fully constrained, no free DOFs)")

    substeps, dangling = parse(a.multi)
    if dangling:
        print(f"NOTE: {dangling} convergence records after the last completed "
              f"substep were dropped (run did not finish cleanly?)")
    for s in substeps:
        print(f"  substep {s['n']}: inc={s['inc']:.4g} "
              f"iters={len(s['rec']) - 1} "
              f"residual {s['rec'][0][0]:.3g} -> {s['rec'][-1][0]:.3g} "
              f"{'CONVERGED' if s['rec'][-1][2] else 'NOT converged'}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(substeps, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
