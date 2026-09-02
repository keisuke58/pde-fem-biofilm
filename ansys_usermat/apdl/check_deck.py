#!/usr/bin/env python3
"""Pre-flight an APDL deck for the failures that produce a plausible wrong answer.

Everything checked here is silent in ANSYS: no error, no warning, a stress
field that looks entirely normal. That is the whole selection criterion --
anything the solver already complains about does not need a script.

    python ansys_usermat/apdl/check_deck.py <deck.dat> [...]

**Time step vs viscous relaxation.** The routine's own guard refuses only
dt/tau > 0.5, where the stress changes sign, because below that the answer
stays qualitatively right and how much accuracy to buy with step size is the
analyst's call. But "qualitatively right" leaves a lot of room: on
t_growth_baseclamped.dat's own material the step AUTOTS is permitted to
coarsen to costs 30% of the von Mises stress. That is the number this thesis
reports.

**State-variable declaration.** usermat_biofilm.f reads the growth driver as
`if (nStatev .ge. 10) ALPHA = ustatev(10)`. A deck that declares fewer leaves
ALPHA at zero, so Fg = I and the solve runs *purely elastic* -- for a thesis
about growth-induced stress, a run with no growth in it, reported as one with
growth. Below 9 the Fv read also runs past the array. Neither says anything.

**Growth driver actually set.** Declaring the slots is not the same as filling
them: a deck with TB,STATE but no TBDATA writing slot 10 has the same silent
no-growth outcome.

**TBDATA value count.** TBDATA takes at most six values per call and drops the
rest without complaint -- a real deck bug in this repository's history.
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
STATEV_FOR_GROWTH = 10     # usermat_biofilm.f: alpha lives in ustatev(10)
STATEV_FOR_FV = 9          # ustatev(1:9) is Fv
STATEV_WITH_MATERIAL = 14  # prop(7)=1 puts C10,C01,D1,eta in ustatev(11:14)
TBDATA_MAX_VALUES = 6      # APDL drops the rest silently


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

    # TB,STATE,mat,ntemp,npts   and   TB,USER,mat,ntemp,npts
    def tb_count(kind):
        m = re.search(rf"^\s*TB\s*,\s*{kind}\s*,\s*\d+\s*,([^!\n]*)", text, re.M | re.I)
        if not m:
            return None
        fields = [f.strip() for f in m.group(1).split(",")]
        nums = [int(f) for f in fields if f.isdigit()]
        return nums[-1] if nums else None

    # every TBDATA call: (start index, how many values)
    tbdata = []
    for m in re.finditer(r"^\s*TBDATA\s*,\s*(\d+)\s*,([^!\n]*)", text, re.M | re.I):
        vals = [v for v in (x.strip() for x in m.group(2).split(",")) if v != ""]
        tbdata.append((int(m.group(1)), len(vals), m.start()))

    # is prop(7) (kStateMat) switched on?
    kstatemat = False
    for start, n, _ in tbdata:
        if start <= 7 < start + n:
            f = [v.strip() for v in re.split(r",", text[_:].split("\n")[0])[2:]]
            try:
                if num(f[7 - start]) > 0.5:
                    kstatemat = True
            except (IndexError, ValueError):
                pass

    return {"c10": c10, "eta": eta, "total": total, "steps": dts,
            "n_state": tb_count("STATE"), "n_prop": tb_count("USER"),
            "tbdata": tbdata, "kstatemat": kstatemat,
            "alpha_written": any(s_ <= STATEV_FOR_GROWTH < s_ + n
                                 for s_, n, _ in tbdata
                                 if s_ <= STATEV_FOR_GROWTH)}


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
        # --- silent-failure checks that have nothing to do with the step ---
        ns = d["n_state"]
        if ns is None:
            print("  no TB,STATE found — the growth driver has nowhere to live")
            worst_ok = False
        else:
            need = STATEV_WITH_MATERIAL if d["kstatemat"] else STATEV_FOR_GROWTH
            if ns < STATEV_FOR_FV:
                print(f"  TB,STATE declares {ns} — below {STATEV_FOR_FV}, the Fv "
                      "read runs past the array")
                worst_ok = False
            elif ns < need:
                print(f"  TB,STATE declares {ns}, needs {need} — ALPHA is never "
                      "read, so Fg = I and this solve is PURELY ELASTIC with no "
                      "warning")
                worst_ok = False
            elif not d["alpha_written"]:
                print(f"  TB,STATE declares {ns} but no TBDATA writes slot "
                      f"{STATEV_FOR_GROWTH} — alpha stays 0, i.e. no growth")
                worst_ok = False

        over = [(st, n) for st, n, _ in d["tbdata"] if n > TBDATA_MAX_VALUES]
        if over:
            for st, n in over:
                print(f"  TBDATA,{st} carries {n} values — APDL takes "
                      f"{TBDATA_MAX_VALUES} and drops the rest silently")
            worst_ok = False

        if d["eta"] <= 0.0:
            print("  eta = 0 (elastic): no relaxation time, no step constraint")
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
        print("\nAt least one deck would produce a plausible wrong answer.")
        print("For a step that is merely too coarse: tighten NSUBST's *third*")
        print("argument (the fewest substeps -- it sets the LARGEST step AUTOTS")
        print("may take), or DELTIM's max.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
