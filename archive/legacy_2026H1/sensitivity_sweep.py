#!/usr/bin/env python3
"""Literature-based sensitivity: sigma_CH/sigma_DH vs nutrient penetration depth Lp.

Physical basis (O2, the usual limiting nutrient in oral biofilm):
  D_eff(O2 in biofilm) ~ 1.5-2.4e-9 m^2/s   (Stewart 2003, J Bacteriol review)
  O2 penetration depth ~ 50-200 um in dense biofilm (de Beer 1994; Stewart 2003)
  biofilm shell thickness L = 0.2 mm (from the tooth mesh)
  -> Lp/L ~ 0.25-1.0 (nutrient-limited); also probe the well-mixed limit.

Lp = sqrt(D_C/GAMMA_C); we set D_C = Lp^2 * GAMMA_C to scan Lp.
Established biofilm phi=1, uniform E (phi-gate=1). CH & DH -> sigma -> ratio.
"""
import sys, subprocess
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import factorized
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "JAXFEM"))
import tooth_pde3d as t3
import gen_tooth_klempt_umat_inp as g

OUT = HERE / "_sens_sweep"; OUT.mkdir(exist_ok=True)
EXTRACT = HERE / "extract_mises.py"
PHI = {"commensal_hobic": np.array([0.942, 0.012, 0.012, 0.011, 0.011]),
       "dysbiotic_hobic": np.array([0.097, 0.119, 0.474, 0.123, 0.093])}
GAMMA_C, K_M = t3.GAMMA_C, t3.K_M

nodes, tets = t3.parse_mesh("p23_biofilm_commensal_static.inp")
Ml, K = t3.assemble(nodes, tets)
outer = np.arange(4 * t3.NV, 5 * t3.NV)
N = len(Ml); M = sp.diags(Ml)
blocks = g.parse_template(g.TMPL)


def alpha_for_Lp(phi_vec, Lp, dt=2e-3, ns=1000):
    D_C = Lp**2 * GAMMA_C
    k_alpha = float(phi_vec @ t3.K_ALPHA)
    A_c = (M + dt * D_C * K).tolil()
    for n in outer:
        A_c.rows[n] = [n]; A_c.data[n] = [1.0]
    sc = factorized(A_c.tocsc())
    c = np.ones(N); alpha = np.zeros(N)
    for _ in range(ns):
        monod = c / (K_M + c)
        rhs = Ml * (c - dt * GAMMA_C * monod); rhs[outer] = 1.0
        c = sc(rhs); np.clip(c, 0, 1, out=c)
        alpha = alpha + dt * (k_alpha * monod)
    return alpha


def sigma(phi_vec, alpha_nodes, tag):
    inp = OUT / ("s_%s.inp" % tag)
    g.write_mode_A(inp, blocks, alpha_nodes, np.ones(N), phi_vec, tag)
    subprocess.run(["abaqus", "job=s_" + tag, "input=" + inp.name, "user=" + str(g.UMAT["A"]),
                    "cpus=1", "interactive", "ask_delete=OFF"],
                   cwd=str(OUT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    r = subprocess.run(["abaqus", "python", str(EXTRACT), str(OUT / ("s_%s.odb" % tag))],
                       cwd=str(OUT), capture_output=True, text=True)
    for ln in r.stdout.splitlines():
        if ln.startswith("MISES_MAX"):
            return float(ln.split()[1]) * 1e3
    return None


L = 0.2
print("Lp[mm]  Lp/L   CH[kPa]  DH[kPa]  ratio   (headline linear-ramp ratio=6.44x)")
for Lp in [0.03, 0.05, 0.10, 0.20, 0.50, 2.0]:
    sigs = {}
    for cond, phi in PHI.items():
        a = alpha_for_Lp(phi, Lp)
        sigs[cond] = sigma(phi, a, "%s_Lp%03d" % (cond[:2], int(Lp * 100)))
    ch, dh = sigs["commensal_hobic"], sigs["dysbiotic_hobic"]
    rr = ch / dh if (ch and dh) else float("nan")
    print("%.2f    %.2f   %7.2f  %7.2f  %5.1fx" % (Lp, Lp / L, ch or -1, dh or -1, rr), flush=True)
