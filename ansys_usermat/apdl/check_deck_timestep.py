#!/usr/bin/env python3
"""Check an APDL deck's time stepping against the viscous relaxation time.

Run this before any solve whose numbers are going to be reported. The
constraint it checks is not enforced anywhere the solver can see: the routine's
own guard refuses only steps past dt/tau = 0.5, where the stress changes sign,
because below that the answer stays qualitatively right and how much accuracy
to buy with step size is the analyst's call. But "qualitatively right" leaves a
lot of room -- measured on t_growth_baseclamped.dat's own material, the step
AUTOTS is permitted to coarsen to costs 30% of the von Mises stress, silently.

That is the number this thesis reports, so it is worth one command.

    python ansys_usermat/apdl/check_deck_timestep.py <deck.dat> [...]

Reads C10 and eta from the deck's TBDATA line and the step range from
TIME/NSUBST/DELTIM, then reports dt/tau and the expected error. Exits non-zero
if any deck can take a step the guard would refuse.
"""
import argparse
import re
import sys
from pathlib import Path

# Measured on the ANSYS core (constrained growth, C10=2e-4, eta=8e-3): the
# von Mises error relative to a fully resolved step, as a function of dt/tau.
# Interpolated between these; they come from ansys_usermat/apdl/ and are
# reproduced by growth_law_verification.ipynb.
_ERROR_CURVE = [(0.001, 0.0), (0.003, 0.2), (0.013, 1.5), (0.062, 7.5),
                (0.125, 15.0), (0.250, 30.0), (0.500, 100.0)]
GUARD_RATIO = 0.5          # DTMAX_RATIO in biofilm_material_v01.f


def expected_error(ratio):
    """Percent von Mises error at this dt/tau, linearly interpolated."""
    if ratio <= _ERROR_CURVE[0][0]:
        return 0.0
    if ratio >= _ERROR_CURVE[-1][0]:
        return float("inf")
    for (r0, e0), (r1, e1) in zip(_ERROR_CURVE, _ERROR_CURVE[1:]):
        if r0 <= ratio <= r1:
            return e0 + (e1 - e0) * (ratio - r0) / (r1 - r0)
    return float("nan")


def parse_deck(path):
    text = path.read_text(errors="replace")

    def num(s):
        return float(s.replace("E", "e").replace("D", "e"))

    tb = re.search(r"^\s*TBDATA\s*,\s*1\s*,(.+)$", text, re.M | re.I)
    if not tb:
        return None
    vals = [v.strip() for v in tb.group(1).split("!")[0].split(",")]
    try:
        c10, eta = num(vals[0]), num(vals[3])
    except (IndexError, ValueError):
        return None

    t = re.search(r"^\s*TIME\s*,\s*([0-9.eEdD+-]+)", text, re.M | re.I)
    total = num(t.group(1)) if t else None

    dts = []
    ns = re.search(r"^\s*NSUBST\s*,\s*([0-9]+)\s*(?:,\s*([0-9]+))?\s*(?:,\s*([0-9]+))?",
                   text, re.M | re.I)
    if ns and total:
        init, mx, mn = (int(g) if g else None for g in ns.groups())
        # NSBMN is the *fewest* substeps, i.e. the LARGEST step AUTOTS may take.
        dts = [("initial", total / init)]
        if mn:
            dts.append(("coarsest (AUTOTS)", total / mn))
        if mx:
            dts.append(("finest (AUTOTS)", total / mx))

    dl = re.search(r"^\s*DELTIM\s*,\s*([0-9.eEdD+-]+)\s*(?:,\s*([0-9.eEdD+-]+))?"
                   r"\s*(?:,\s*([0-9.eEdD+-]+))?", text, re.M | re.I)
    if dl:
        dts.append(("DELTIM initial", num(dl.group(1))))
        if dl.group(3):
            dts.append(("DELTIM max", num(dl.group(3))))

    return {"c10": c10, "eta": eta, "total": total, "steps": dts}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("decks", nargs="+", type=Path)
    args = ap.parse_args(argv)

    worst_ok = True
    for deck in args.decks:
        d = parse_deck(deck)
        print(f"\n{deck.name}")
        if d is None:
            print("  no TBDATA line found — skipped")
            continue
        if d["eta"] <= 0.0:
            print(f"  eta = 0 (elastic): no relaxation time, no step constraint")
            continue
        if d["c10"] <= 0.0:
            print("  C10 = 0 — cannot form a relaxation time")
            continue
        tau = d["eta"] / (2.0 * d["c10"])
        print(f"  C10 = {d['c10']:g}, eta = {d['eta']:g}  ->  tau = {tau:g} s")
        if not d["steps"]:
            print("  no TIME/NSUBST/DELTIM found — check the step manually")
            continue
        for label, dt in d["steps"]:
            ratio = dt / tau
            err = expected_error(ratio)
            if ratio > GUARD_RATIO:
                mark, worst_ok = "REFUSED by the routine's guard", False
            elif err >= 5.0:
                mark = f"~{err:.0f}% von Mises error — too coarse to report"
                worst_ok = False
            elif err >= 1.0:
                mark = f"~{err:.1f}% error"
            else:
                mark = "ok"
            print(f"    {label:20s} dt = {dt:9.4g}  dt/tau = {ratio:7.4f}   {mark}")

    if not worst_ok:
        print("\nAt least one permitted step is too coarse for reported numbers.")
        print("Tighten NSUBST's *third* argument (the minimum substep count --")
        print("it sets the LARGEST step AUTOTS may take), or use DELTIM's max.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
