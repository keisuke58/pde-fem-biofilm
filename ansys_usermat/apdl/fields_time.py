"""fields_time.py -- collect fields_1..4.csv (post_fields_time.mac) into
results/2026-10-cal5/fields_<run>.json (A-15): {"run", "t": [4 set times],
"elem", "cx", "cy", "cz", "phi": [[per time]], "nut", "alpha", "seqv"},
7 significant digits. phi = SVAR 41 (field phi at the Gauss points, element
mean), nut = SVAR 51, alpha = SVAR 84 (alpha - 1), seqv in the deck's units.

    python ansys_usermat/apdl/fields_time.py RUN DIR_WITH_CSV [N]   (N = 40 with post_fields_dense.mac
    -> fields_<run>_dense.json)
"""
import csv
import json
import sys
from pathlib import Path


def g7(x):
    return float(f"{x:.7g}")


def main(run, d, n=4):
    out = {"run": run, "t": []}
    for k in range(1, n + 1):
        rows = list(csv.DictReader(open(Path(d) / f"fields_{k}.csv")))
        rows = [{c.strip(): v for c, v in r.items()} for r in rows]
        rows = [r for r in rows if float(r["cx"]) != 0 or float(r["cy"]) != 0 or float(r["cz"]) != 0
                or float(r["phi"]) != 0]
        out["t"].append(g7(float(rows[0]["t"])))
        if k == 1:
            out["elem"] = [int(float(r["elem"])) for r in rows]
            for c in ("cx", "cy", "cz"):
                out[c] = [g7(float(r[c])) for r in rows]
        else:
            assert [int(float(r["elem"])) for r in rows] == out["elem"]
        for c in ("phi", "nut", "alpha", "seqv"):
            out.setdefault(c, []).append([g7(float(r[c])) for r in rows])
    p = Path(__file__).parent / "results" / "2026-10-cal5" / (f"fields_{run}.json" if n == 4 else f"fields_{run}_dense.json")
    p.write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote {p}: t = {out['t']}, {len(out['elem'])} elements")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
