"""make_wp3_decks.py -- the ANSYS runs of the Keio WP3 study (NEM against
Galerkin, RUN_WP3_IKMHIWI03.md): writes the decks with make_wired_deck.py /
refine_deck.py and one run list per block for run_chain.ps1.

    python ansys_usermat/apdl/make_wp3_decks.py [--block A|B|both] [--check]

Block A, seed growth without the front term (F:\\biofilm_upf_wired, the week
chain's executable): the time steps not run yet next to ds{8,16,24}_beta002*.
Base decks ds8/ds24_beta002_dt4 (beta = 0.02 mm^2/T*, k_alpha = 1e-3, E = 10 Pa,
T* = 1.1); only deltim changes.

Block B, Klempt 2024 test cases 4.2 high and 4.1 with the front term
(F:\\biofilm_upf_front, FRONT_TERM_FIX.md "Front build for WP3"): base
ds_fig7h (Table 2 converted to the 2 mm cube, consumption g phi, T* = 1),
refined to 16^3 / 24^3 with refine_deck.py. 4.1 replaces the components:
NUTRIENT1 = the corner element at (max x, max y, max z), BIOFILM1 = the
elements whose centroid lies within L/4 of the centre (the rule of
abaqus_composition/partner_elem_sets.py --case 41). dt 0.01 / 0.005 / 0.0025
(front Courant number v dt / h <= 0.6 with v <= 5 mm/T*).

Run lists: <workdir>\\_wp3_<block>_runs.txt, one "deck::minutes" per line, in
run order (coarse meshes first). --check prints them and writes nothing.
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WIRED = Path(r"F:\biofilm_upf_wired")
FRONT = Path(r"F:\biofilm_upf_front")


def num(x):
    return f"{x:g}".replace(".", "")


def block_a():
    """(name, base, dt, minutes); base in WIRED."""
    return [("wp3_n8_dt00125", "ds8_beta002_dt4", 0.0125, 20),
            ("wp3_n24_dt01", "ds24_beta002_dt4", 0.1, 90),
            ("wp3_n24_dt005", "ds24_beta002_dt4", 0.05, 150),
            ("wp3_n24_dt00125", "ds24_beta002_dt4", 0.0125, 450)]


# timeouts: about 1.6 x the estimate of RUN_WP3_IKMHIWI03.md
MINUTES_B = {8: {0.01: 20, 0.005: 30, 0.0025: 60},
             16: {0.01: 60, 0.005: 120, 0.0025: 240},
             24: {0.01: 420, 0.005: 840, 0.0025: 1700}}


def block_b():
    """(name, case, mesh, dt, minutes); coarse meshes first."""
    out = []
    for n in (8, 16, 24):
        for case in ("42h", "41"):
            for dt in (0.01, 0.005, 0.0025):
                out.append((f"wp3_{case}_n{n}_dt{num(dt)}", case, n, dt, MINUTES_B[n][dt]))
    return out


def run(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(map(str, cmd))}\n{r.stdout}{r.stderr}")
    return r.stdout


def case41_sets(deck: Path):
    spec = importlib.util.spec_from_file_location(
        "pes", REPO / "abaqus_composition" / "partner_elem_sets.py")
    pes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pes)
    nodes, elems, _ = pes.read(deck)
    cen = {e: [sum(nodes[k][i] for k in c) / 8 for i in range(3)] for e, c in elems.items()}
    lo = [min(v[i] for v in nodes.values()) for i in range(3)]
    hi = [max(v[i] for v in nodes.values()) for i in range(3)]
    L = hi[0] - lo[0]
    mid = [(p + q) / 2 for p, q in zip(lo, hi)]
    nut = [max(cen, key=lambda e: sum(cen[e]))]
    seed = sorted(e for e, c in cen.items()
                  if sum((c[i] - mid[i]) ** 2 for i in range(3)) <= (0.25 * L) ** 2 + 1e-12)
    return nut, seed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--block", choices=("A", "B", "both"), default="both")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    py = sys.executable
    mk = HERE / "make_wired_deck.py"

    if a.block in ("A", "both"):
        lines = [f"{name}::{mins}" for name, _, _, mins in block_a()]
        if not a.check:
            for name, base, dt, _ in block_a():
                b = (WIRED / f"{base}.dat").read_text(errors="replace")
                pe = re.search(r"^\*USE,post_elem_stress\.mac,(\d+)", b, re.M)   # as make_week_decks.py
                run([py, mk, WIRED / f"{base}.dat", WIRED / f"{name}.dat", "--deltim", f"{dt:g}",
                     "--post", "both", "--post-elem", pe.group(1) if pe else "220"])
            (WIRED / "_wp3_A_runs.txt").write_text("\n".join(lines) + "\n")
        print("# block A (" + str(WIRED) + ")\n" + "\n".join(lines))

    if a.block in ("B", "both"):
        rows = block_b()
        lines = [f"{name}::{mins}" for name, _, _, _, mins in rows]
        if not a.check:
            base = {8: FRONT / "ds_fig7h.dat"}
            for n in (16, 24):
                base[n] = FRONT / f"ds_fig7h_n{n}.dat"
                run([py, HERE / "refine_deck.py", base[8], base[n], "--n", str(n)])
            sets = {}
            for n in (8, 16, 24):
                nut, seed = case41_sets(base[n])
                sets[n] = ["--cmblock", "NUTRIENT1=" + ",".join(map(str, nut)),
                           "--cmblock", "BIOFILM1=" + ",".join(map(str, seed))]
                print(f"# 4.1 on {n}^3: nutrient element {nut}, seed {len(seed)} elements")
            for name, case, n, dt, _ in rows:
                cmd = [py, mk, base[n], FRONT / f"{name}.dat", "--deltim", f"{dt:g}", "--time", "1.0",
                       "--post", "all"]
                if case == "41":
                    cmd += sets[n]
                run(cmd)
            (FRONT / "_wp3_B_runs.txt").write_text("\n".join(lines) + "\n")
        print("# block B (" + str(FRONT) + ")\n" + "\n".join(lines))
        hours = sum(m for *_, m in rows) / 60
        print(f"# block B timeouts add up to {hours:.0f} h", file=sys.stderr)


if __name__ == "__main__":
    main()
