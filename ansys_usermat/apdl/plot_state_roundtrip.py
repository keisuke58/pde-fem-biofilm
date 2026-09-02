#!/usr/bin/env python3
"""What comes back out of ANSYS's TB,STATE slots.

ANSYS has no INISTATE for a USERMAT, so the growth variable cannot be given as
a field. It is carried in a state slot instead: nine slots hold Fv row-major
and the tenth holds alpha, written per material with TBDATA. That makes the
state array the single point of failure for the whole growth mechanism, and it
fails silently in both directions:

  * TBDATA takes at most six values per call, and drops the rest without an
    error. The decks therefore split the write across three calls, and a lost
    tail would leave alpha at zero -- an entirely successful, entirely elastic
    run that reports no problem.
  * TB,STATE with too few slots has the same effect. check_deck.py catches
    both before a run; this script confirms the outcome after one.

Reading the slots back is the only direct evidence the write survived. This
script does that for the three v222 runs that printed SVAR.

Two things it does NOT establish, stated on the figure so the reader is not
left to assume otherwise. Every deck sets eta = 0, so Fv has nothing to
evolve: the round trip is exercised, the viscous update is not. And the
cylinder deck lists SVAR under ESEL,S,MAT,,2, so it evidences the growth layer
only -- the substrate's alpha = 0 is not read back anywhere.

That ESEL is not incidental. An earlier version of the deck had both volumes
land on the default material, where ESEL,S,MAT,,2 selected zero elements; the
11040 it selects here is itself the check that the two-material assignment
took.

    python ansys_usermat/apdl/plot_state_roundtrip.py -o assets/v222_state_roundtrip.png
"""
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

SLOTS = 10
LABELS = [r"$F_v$" + f"[{i//3+1},{i%3+1}]" for i in range(9)] + [r"$\alpha$"]
IDENT = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]

_BLOCK = re.compile(r"ELEMENT=\s*(\d+)\s+SOLID185\s*\n\s*NODE\s+(\d+)\s*\n"
                    r"((?:\s+\d+\s+[-+0-9.E]+\s*\n)+)")
_ROW = re.compile(r"\s+(\d+)\s+([-+0-9.E]+)\s*\n")


def read_slots(path):
    """Returns {slot: (n_samples, {value: count})} for every SVAR block found.

    Values are compared as printed. MAPDL lists five significant digits, so
    "matches" here means matches to the precision the solver reports -- not
    bit equality, which the listing cannot show.
    """
    txt = Path(path).read_text(errors="replace")
    per = {}
    for m in _BLOCK.finditer(txt):
        slot = int(m.group(2))
        vals = [float(v) for _, v in _ROW.findall(m.group(3))]
        n, c = per.get(slot, (0, Counter()))
        c.update(round(v, 12) for v in vals)
        per[slot] = (n + len(vals), c)
    return per


def plot(runs, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10.2, 3.2))
    ok, absent = "#1e8449", "0.88"

    for r, (name, expect_alpha, per) in enumerate(runs):
        y = len(runs) - 1 - r
        for s in range(1, SLOTS + 1):
            want = IDENT[s - 1] if s <= 9 else expect_alpha
            if s not in per:
                ax.add_patch(plt.Rectangle((s - 1.47, y - 0.38), 0.94, 0.76,
                                           fc=absent, ec="white", lw=1.4))
                ax.text(s - 1, y, "not\nlisted", ha="center", va="center",
                        fontsize=6.5, color="0.45")
                continue
            n, c = per[s]
            good = len(c) == 1 and abs(next(iter(c)) - want) <= 1e-12
            ax.add_patch(plt.Rectangle((s - 1.47, y - 0.38), 0.94, 0.76,
                                       fc=ok if good else "#c0392b",
                                       ec="white", lw=1.4))
            got = next(iter(c)) if len(c) == 1 else None
            txt = f"{got:g}" if got is not None else "MIXED"
            ax.text(s - 1, y + 0.09, txt, ha="center", va="center",
                    fontsize=8.5, color="white")
            ax.text(s - 1, y - 0.19, f"{n:,} samples", ha="center", va="center",
                    fontsize=6.2, color="white")

    ax.set_xlim(-0.55, SLOTS - 0.45)
    ax.set_ylim(-0.6, len(runs) - 0.4)
    ax.set_xticks(range(SLOTS))
    ax.set_xticklabels(LABELS, fontsize=8.5)
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels([r[0] for r in reversed(runs)], fontsize=8.5)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title("TB,STATE round trip: what was written, read back out of the "
                 "solve", fontsize=11)
    fig.text(0.045, 0.145,
             "green = every sample equals what the deck wrote, to the five "
             "digits MAPDL prints.\n"
             r"$\eta = 0$ in all three decks, so $F_v$ has nothing to evolve: "
             "the round trip is exercised, the viscous update is not.\n"
             "The cylinder lists SVAR under ESEL,S,MAT,,2, so it evidences the "
             r"growth layer only -- the substrate's $\alpha = 0$ is not read "
             "back anywhere.",
             fontsize=7.8, color="0.35", va="top", ha="left")
    fig.tight_layout(rect=(0, 0.20, 1, 1))
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


def main(argv=None):
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("assets/v222_state_roundtrip.png"))
    a = ap.parse_args(argv)

    spec = [("constrained, 1 element", 0.05, "growth_result.txt"),
            ("free, 1 element", 0.05, "growth_free_result.txt"),
            ("curved shell, growth layer\n(11040 of 12240 elements)", 0.01,
             "growth_cylinder_wrapper_result.txt")]
    runs, bad = [], 0
    for name, alpha, fn in spec:
        per = read_slots(here / fn)
        runs.append((name, alpha, per))
        for s, (n, c) in sorted(per.items()):
            want = IDENT[s - 1] if s <= 9 else alpha
            got = next(iter(c)) if len(c) == 1 else None
            flag = "" if got is not None and abs(got - want) <= 1e-12 else "  <-- MISMATCH"
            bad += bool(flag)
            print(f"{fn:42s} slot {s:2d}: {n:6d} samples, "
                  f"wrote {want:g}, read {got if got is not None else dict(c)}{flag}")
    if bad:
        print(f"{bad} slot(s) did not read back what the deck wrote")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    plot(runs, a.out)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
