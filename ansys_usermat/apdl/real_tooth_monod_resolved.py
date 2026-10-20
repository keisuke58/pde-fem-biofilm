"""real_tooth_monod_resolved.py -- the converged stress of the end states of the Monod runs on the real molar
(real_tooth_model.py monod, monod_both, monod_r01), with the film submodel of real_tooth_submodel.py at h = h_f / 2
(rigid tooth, crop box, void phi <= 0.05 dropped, voxel wall), as real_tooth_physics_resolved.py does for the other
runs. The Monod films are thicker than 0.5 mm, so the resolution across the film is better than the six elements
checked there.

    python ansys_usermat/apdl/real_tooth_monod_resolved.py [tags]
        -> merged into ansys_usermat/apdl/results/2026-10-3d/real_tooth_physics_resolved.json
"""
import json
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import real_tooth_submodel as SM  # noqa: E402

OUT = HERE / "results" / "2026-10-3d" / "real_tooth_physics_resolved.json"
tags = sys.argv[1:] or ["monod_r01", "monod"]
for tag in tags:
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    if tag in res:
        continue
    F = SM.load_field(tag)
    r = SM.solve_submodel(F["phi"], F["am1"], F["h_f"] / 2, tooth=F["tooth"], sdf=F["sdf"], h_f=F["h_f"],
                          void="drop", void_thr=0.05, tooth_mode="voxel", verbose=True)
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res[tag] = dict(measures=r["measures"], seconds=r["seconds"], h_mm=F["h_f"] / 2, source=F["source"])
    OUT.write_text(json.dumps(res, indent=1, default=float))
    m = r["measures"]
    print(tag, {k: m[k] for k in ("film", "band", "core", "band/core", "p mean") if k in m}, flush=True)
