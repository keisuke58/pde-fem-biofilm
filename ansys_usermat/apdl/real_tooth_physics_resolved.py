"""real_tooth_physics_resolved.py -- the converged stress of the final states of the real-molar runs.

real_tooth_model.py solved the stress on the 32^3 full-box mesh (h = 0.5 mm, the 0.5 mm film one element thick), which
REAL_TOOTH_CONVERGENCE.md shows is not converged. Here the final state of every run (96^3 field) is solved again with
the film submodel of real_tooth_submodel.py (rigid tooth, crop box, void phi <= 0.05 dropped, voxel wall) at
h = h_f / 2 = 0.083 mm, six elements across the film, where REAL_TOOTH_SUBMODEL.md finds the film mean within 1 % of the
finest meshes. The gate runs were driven by the 32^3 pressure during the run; only their end state is re-solved here.

    python ansys_usermat/apdl/real_tooth_physics_resolved.py
        -> ansys_usermat/apdl/results/2026-10-3d/real_tooth_physics_resolved.json
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
RUNS = ["base", "nutrient", "gate", "nutrient_gate", "nutrient_t2"]
res = json.loads(OUT.read_text()) if OUT.exists() else {}
for tag in RUNS:
    if tag in res:
        continue
    F = SM.load_field(tag)
    r = SM.solve_submodel(F["phi"], F["am1"], F["h_f"] / 2, tooth=F["tooth"], sdf=F["sdf"], h_f=F["h_f"],
                          void="drop", void_thr=0.05, tooth_mode="voxel", verbose=False)
    res[tag] = dict(measures=r["measures"], seconds=r["seconds"], h_mm=F["h_f"] / 2, source=F["source"])
    OUT.write_text(json.dumps(res, indent=1, default=float))
    m = r["measures"]
    print(tag, {k: m[k] for k in ("film", "band", "core", "band/core", "p mean") if k in m}, flush=True)
