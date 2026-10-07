#!/usr/bin/env python3
"""Is Klempt 2024's factor of ten a resolution artefact?

klempt2024_sensitivity.py cleared the seed and the nutrient by measurement,
and reading the paper cleared the constants and showed that the variants which
fit best are not the paper's. What is left, named in
KLEMPT2024_REPRODUCTION.md and never tested, is the scheme: an explicit upwind
finite-difference solve here against the paper's implicit Galerkin FEM.

There is a specific reason to suspect resolution rather than the solver. The
growth term is proportional to ||grad phi||, so whatever under-resolves the
interface under-states the growth directly and by the same factor. On a 1 um
grid with beta = 2 um^2/T* and a growth rate of order r c/(k+c) ~ 50 /T*, the
interface the phase field wants is around sqrt(beta / rate) ~ 0.2 um -- five
times finer than the grid it is being resolved on.

So this refines the grid and watches, and it runs THE PAPER'S OWN VARIANT
(growth = the signed dot product of Eq. 34, consumption = g phi of Eq. 35 and
Eq. 24), not the combination that happens to fit. Anything else would be
measuring how a variant the paper does not use responds to refinement.

If phi(T*) climbs toward the paper's curve as h falls, the discrepancy is
discretisation and this repository can close it. If it converges somewhere
short, the discrepancy is not resolution and what remains is the deforming
domain, the clip on phi, or the equations not being what produced the figures
-- and that last one is a question for the authors.

N is a module-level constant in klempt2024_quantitative (the grid, the
weights and the Laplacian are all built from it at import), so each resolution
is a fresh execution of that module's source with N substituted, rather than a
parameter.

Cost grows steeply: the nutrient solve is sparse-direct on N^3 unknowns, so
N = 41 is roughly fifteen times N = 21. Start small.

    python JAXFEM/klempt2024_resolution.py --grids 21 31
"""
import argparse
import json
import sys
import time
import types
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_SRC = _HERE / "klempt2024_quantitative.py"

PAPER_T = [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00]


def load_at(n):
    """klempt2024_quantitative with N = n, as its own module instance."""
    src = _SRC.read_text()
    old = "\nN = 21\n"
    if old not in src:
        raise SystemExit("could not find 'N = 21' in klempt2024_quantitative.py "
                         "-- the grid constant moved, fix this loader")
    src = src.replace(old, f"\nN = {n}\n", 1)
    mod = types.ModuleType(f"k24_N{n}")
    mod.__file__ = str(_SRC)
    sys.path.insert(0, str(_HERE))
    exec(compile(src, str(_SRC), "exec"), mod.__dict__)
    return mod


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--grids", type=int, nargs="+", default=[21, 31])
    ap.add_argument("--case", default="fig7_high")
    ap.add_argument("--t-end", type=float, default=0.25)
    ap.add_argument("-o", "--out", type=Path,
                    default=Path("JAXFEM/klempt2024_results/resolution.json"))
    a = ap.parse_args(argv)

    ref = load_at(21)
    paper = ref.PAPER[a.case]
    ts = [t for t in PAPER_T if t <= a.t_end + 1e-12]
    print(f"case {a.case}, the paper's own variant "
          f"(growth=printed, consumption=printed), to T* = {a.t_end}")
    head = "  ".join(f"{t:5.2f}" for t in ts)
    print(f"{'':22s}{head}")
    print(f"{'PAPER':22s}" +
          "  ".join(f"{v:5.3f}" for v, t in zip(paper['phi'], paper['t']) if t in ts))

    out, t0 = {}, time.time()
    for n in a.grids:
        m = load_at(n)
        rec = m.run(a.case, "printed", growth="printed", t_end=a.t_end)
        phi = np.interp(ts, rec["t"], rec["phi"])
        out[str(n)] = {"h": m.H, "t": ts, "phi": [round(float(v), 4) for v in phi]}
        print(f"N={n:<3d} h={m.H:.3f} um     " +
              "  ".join(f"{v:5.3f}" for v in phi), flush=True)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(
        {"case": a.case, "variant": "paper (growth=printed, consumption=printed)",
         "t_end": a.t_end, "paper": {"t": paper["t"], "phi": paper["phi"]},
         "grids": out}, indent=1))
    print(f"\nwrote {a.out}  [{time.time()-t0:.0f} s]")
    if len(out) >= 2:
        ks = sorted(out, key=int)
        lo, hi = out[ks[0]]["phi"][-1], out[ks[-1]]["phi"][-1]
        tgt = np.interp(ts[-1], paper["t"], paper["phi"])
        print(f"at T* = {ts[-1]}: h={out[ks[0]]['h']:.2f} gives {lo:.3f}, "
              f"h={out[ks[-1]]['h']:.2f} gives {hi:.3f}, the paper {tgt:.3f}")
        if hi > lo * 1.2:
            print("  -> refinement moves it toward the paper: discretisation is "
                  "part of the gap")
        else:
            print("  -> refinement does not move it: the gap is not resolution")
    return 0


if __name__ == "__main__":
    sys.exit(main())
