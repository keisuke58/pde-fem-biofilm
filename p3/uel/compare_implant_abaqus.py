"""compare_implant_abaqus.py -- Abaqus (Gauss-point file of the last increment) against the Python emulation for one
implant-sector mesh.

    python compare_implant_abaqus.py <tag> [gp_file]       tag as in results/implant_<tag>.json, e.g. 4x12x16_b1
Prints and writes results/implant_<tag>_abaqus_check.json: largest difference of alpha - 1 and vM over the elements (Gauss-point means), relative to the
largest value of the field (elements with alpha - 1 = 0 make a pointwise ratio undefined), and the seed-layer means of both.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
tag = sys.argv[1]
gp = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "abaqus_in" / f"gp_{tag}.txt"
blk = []
for line in open(gp):
    if line.startswith("# step"):
        blk = []
    elif not line.startswith("#") and line.strip():
        blk.append(line.split())
g = np.array(blk, float)
el = g[:, 0].astype(int) - 1
ne = el.max() + 1
mean = lambda c: np.array([g[el == e, c].mean() for e in range(ne)])  # noqa: E731
a, v = mean(7), mean(9)
P = json.load(open(HERE / "results" / f"implant_{tag}.json"))
pa, pv = np.array(P["elements"]["alpha_minus_1"]), np.array(P["elements"]["vM_over_E"])
seed = np.array(P["elements"]["layer"]) == 0
out = dict(tag=tag, elements=int(ne), max_rel_alpha=float(np.abs(a - pa).max() / np.abs(pa).max()), max_rel_vM=float(np.abs(v - pv).max() / np.abs(pv).max()),
           seed_alpha_abaqus=float(a[seed].mean()), seed_alpha_python=float(pa[seed].mean()),
           seed_vM_abaqus=float(v[seed].mean()), seed_vM_python=float(pv[seed].mean()),
           max_vM_abaqus=float(v.max()), max_vM_python=float(pv.max()))
(HERE / "results" / f"implant_{tag}_abaqus_check.json").write_text(json.dumps(out, indent=1))
print(out)
