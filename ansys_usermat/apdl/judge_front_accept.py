"""judge_front_accept.py -- acceptance of the front-term build
(FRONT_TERM_FIX.md, "Front build for WP3") from three runs of ds_fig7h (8^3,
T* = 1) in F:\\biofilm_upf_front:

    fr_fig7h_dt0005    front term, dt 0.005
    fr_fig7h_dt00025   front term, dt 0.0025
    fr_fig7h_nofront   MAX_GROWTH11 = 0, dt 0.005

    python ansys_usermat/apdl/judge_front_accept.py [--workdir F:\\biofilm_upf_front] [--report file]

Reads the FRONTPHI lines (time, mean, min, max of phi per substep) from
out_<run>.txt and alpha - 1 per element (all_stress_<run>.csv, SVAR 84, a
time integral of phi) for the direction check. Two levels:

  strict (FRONT_TERM_FIX.md): 0 <= phi <= 1 within 1e-6; mean phi at dt and
          dt/2 within 1 %
  gate   (whether the WP3 runs of block B start): all three runs reach
          T* = 1 with finite values; -0.05 <= phi <= 1.05; mean phi at dt
          and dt/2 within 5 %; direction: the extra alpha - 1 from the front
          term on the nutrient side of the seed is positive and at least
          twice that on the far side.

Exit code 0 when the gate passes, 1 otherwise. The strict result is reported,
it does not stop the runs.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

RUNS = ("fr_fig7h_dt0005", "fr_fig7h_dt00025", "fr_fig7h_nofront")


def front_phi(path: Path):
    rows = {}
    for line in open(path, errors="replace"):
        if line.startswith("FRONTPHI"):
            v = [float(x) for x in line.split()[1:5]]
            rows[v[0]] = v
    return [rows[t] for t in sorted(rows)]


def alpha_by_elem(path: Path):
    out = {}
    with open(path) as f:
        r = csv.reader(f)
        head = [h.strip() for h in next(r)]
        ia, iy = head.index("alpha"), head.index("cy")
        for row in r:
            if row:
                out[int(float(row[0]))] = (float(row[ia]), float(row[iy]))
    return out


def component_ids(deck: Path, name: str):
    lines = deck.read_text(errors="replace").splitlines()
    for i, l in enumerate(lines):
        if re.match(rf"CMBLOCK,{name}\s*,ELEM", l, re.I):
            n, ids, j = int(l.split(",")[3]), [], i + 2
            while len(ids) < n:
                ids += [int(x) for x in re.findall(r"-?\d+", lines[j])]
                j += 1
            full = []
            for v in ids:
                full += list(range(full[-1] + 1, -v + 1)) if v < 0 else [v]
            return full
    raise SystemExit(f"{name} not found in {deck}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=r"F:\biofilm_upf_front")
    ap.add_argument("--report")
    a = ap.parse_args()
    W = Path(a.workdir)
    msg, gate, strict = [], True, True

    def say(s=""):
        msg.append(s)
        print(s)

    say("Front-term build, acceptance on ds_fig7h (8^3, T* = 1)")
    say()
    hist = {}
    for run in RUNS:
        out, st = W / f"out_{run}.txt", W / f"all_stress_{run}.csv"
        h = front_phi(out) if out.exists() else []
        ok = bool(h) and h[-1][0] >= 0.99 and st.exists() and all(math.isfinite(x) for row in h for x in row)
        hist[run] = h
        say(f"{run}: {'complete' if ok else 'NOT complete'}"
            + (f", last T* {h[-1][0]:.4f}, {len(h)} substeps" if h else ", no FRONTPHI lines"))
        gate &= ok
    if not gate:
        say()
        say("GATE: FAIL (a run did not finish)")
        return finish(a, msg, 1)

    f1, f2, f0 = (hist[r] for r in RUNS)
    say()
    say("T*      mean phi: dt 0.005   dt 0.0025   no front    min phi (dt 0.005)  max phi (dt 0.005)")
    for t in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0):
        pick = [min(h, key=lambda r: abs(r[0] - t)) for h in (f1, f2, f0)]
        say(f"{t:4.2f}   {pick[0][1]:12.5f} {pick[1][1]:11.5f} {pick[2][1]:11.5f}   {pick[0][2]:18.3e} {pick[0][3]:18.4f}")

    lo = min(r[2] for h in (f1, f2) for r in h)
    hi = max(r[3] for h in (f1, f2) for r in h)
    b_strict = lo >= -1e-6 and hi <= 1 + 1e-6
    b_gate = lo >= -0.05 and hi <= 1.05
    say()
    say(f"bounds: min phi {lo:.3e}, max phi {hi:.5f}  strict {'pass' if b_strict else 'fail'}, gate {'pass' if b_gate else 'fail'}")

    m1, m2 = f1[-1][1], f2[-1][1]
    rel = abs(m1 - m2) / max(abs(m2), 1e-30)
    worst = max(abs(min(f1, key=lambda r: abs(r[0] - t))[1] - min(f2, key=lambda r: abs(r[0] - t))[1])
                / max(abs(min(f2, key=lambda r: abs(r[0] - t))[1]), 1e-30) for t in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0))
    d_strict, d_gate = worst <= 0.01, worst <= 0.05
    say(f"time step: mean phi at T* = 1 {m1:.5f} / {m2:.5f} ({100 * rel:.2f} %), largest difference over T* {100 * worst:.2f} %"
        f"  strict {'pass' if d_strict else 'fail'}, gate {'pass' if d_gate else 'fail'}")

    deck = W / f"{RUNS[0]}.dat"
    a1, a0 = alpha_by_elem(W / f"all_stress_{RUNS[0]}.csv"), alpha_by_elem(W / f"all_stress_{RUNS[2]}.csv")
    seed, nut = component_ids(deck, "BIOFILM1"), component_ids(deck, "NUTRIENT1")
    ys = sum(a1[e][1] for e in seed) / len(seed)
    yn = sum(a1[e][1] for e in nut) / len(nut)
    sgn = 1.0 if yn > ys else -1.0
    near = [e for e in a1 if sgn * (a1[e][1] - ys) > 0.2 and e not in seed]
    far = [e for e in a1 if sgn * (a1[e][1] - ys) < -0.2 and e not in seed]
    gn = sum(a1[e][0] - a0[e][0] for e in near) / max(len(near), 1)
    gf = sum(a1[e][0] - a0[e][0] for e in far) / max(len(far), 1)
    dir_gate = gn > 0 and gn >= 2 * max(gf, 0.0)
    say(f"direction: seed at y = {ys:.3f}, nutrient at y = {yn:.3f}; extra alpha - 1 from the front term,"
        f" nutrient side {gn:.3e} ({len(near)} elements), far side {gf:.3e} ({len(far)} elements)"
        f"  gate {'pass' if dir_gate else 'fail'}")

    strict = b_strict and d_strict
    gate = b_gate and d_gate and dir_gate
    say()
    say(f"STRICT: {'PASS' if strict else 'FAIL'}")
    say(f"GATE: {'PASS' if gate else 'FAIL'}")
    return finish(a, msg, 0 if gate else 1)


def finish(a, msg, rc):
    if a.report:
        Path(a.report).parent.mkdir(parents=True, exist_ok=True)
        Path(a.report).write_text("\n".join(msg) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
