#!/usr/bin/env python3
"""Run the UMAT with alpha(node) from the 3D FE PDE (no depth ramp / no mapping).
CH & DH -> sigma -> ratio, compared to the headline (linear-ramp) 6.44x."""
import sys, subprocess, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "JAXFEM"))
import tooth_pde3d as t3
import gen_tooth_klempt_umat_inp as g

OUT = HERE / "_pde3d_umat"; OUT.mkdir(exist_ok=True)
PHI = {"commensal_hobic": np.array([0.942, 0.012, 0.012, 0.011, 0.011]),
       "dysbiotic_hobic": np.array([0.097, 0.119, 0.474, 0.123, 0.093])}
MAP = {"commensal_hobic": 13.72, "dysbiotic_hobic": 2.13}

nodes, tets = t3.parse_mesh("p23_biofilm_commensal_static.inp")
Ml, K = t3.assemble(nodes, tets)
outer = np.arange(4 * t3.NV, 5 * t3.NV)
blocks = g.parse_template(g.TMPL)
EXTRACT = HERE / "extract_mises.py"

sig = {}
for cond, phi in PHI.items():
    alpha_nodes, phi_field = t3.solve(Ml, K, phi, outer)   # per-node, 3D PDE
    # phi_local for Klempt Eq.20 gating: use the 3D phi field directly (no max*depth)
    phi_nodes = np.clip(phi_field, 0.0, 1.0)
    inp = OUT / ("pde3d_%s.inp" % cond)
    g.write_mode_A(inp, blocks, alpha_nodes, phi_nodes, phi, cond)
    job = "pde3d_" + cond
    subprocess.run(["abaqus", "job=" + job, "input=" + inp.name, "user=" + str(g.UMAT["A"]),
                    "cpus=1", "interactive", "ask_delete=OFF"],
                   cwd=str(OUT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    r = subprocess.run(["abaqus", "python", str(EXTRACT), str(OUT / (job + ".odb"))],
                       cwd=str(OUT), capture_output=True, text=True)
    s = None
    for ln in r.stdout.splitlines():
        if ln.startswith("MISES_MAX"):
            s = float(ln.split()[1]) * 1e3   # kPa
    sig[cond] = s
    print("[%s] 3D-PDE: alpha_max=%.3f mean=%.3f  sigma=%.2f kPa  (headline linear-ramp=%.2f)"
          % (cond, alpha_nodes.max(), alpha_nodes.mean(), s if s else -1, MAP[cond]), flush=True)

if sig["commensal_hobic"] and sig["dysbiotic_hobic"]:
    r3d = sig["commensal_hobic"] / sig["dysbiotic_hobic"]
    print("\n  sigma_CH/sigma_DH :  3D-PDE = %.2fx   vs   headline(linear-ramp) = 6.44x" % r3d)
