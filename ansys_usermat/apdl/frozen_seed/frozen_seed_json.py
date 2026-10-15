"""frozen_seed_json.py -- the ANSYS side of request A-16 in the form of
p3/uel/results/seed_frozen_<n>.json: per element centroid, mean nodal phi,
alpha - 1 and von Mises (Gauss-point means, ETABLE), p = (SX + SY + SZ)/3,
the seed means (mean nodal phi > 0.5) and the largest von Mises; plus two
checks: the Gauss-point phi against the trilinear interpolation of the nodal
phi, and the stiffness factor recorded at the Gauss points against phi^2 + f.

    python frozen_seed_json.py WORKDIR N
"""
import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "p3" / "uel"))
from make_seed_inp import seed_phi  # noqa: E402

E, NU, FV, KA, T, NSTEP = 10.0, 0.49, 1e-3, 1e-3, 1.1, 11


def main(wd, n):
    wd, n = Path(wd), int(n)
    h = 1.0 / n
    rows = [{k.strip(): float(v) for k, v in r.items()} for r in csv.DictReader(open(wd / "frozen_elems.csv"))]
    cen = np.array([[r["cx"], r["cy"], r["cz"]] for r in rows])
    corners = np.array(list(itertools.product((-0.5, 0.5), repeat=3))) * h
    pe = np.array([seed_phi(c + corners).mean() for c in cen])
    xi = np.array([r["alpha"] for r in rows])
    vm = np.array([r["seqv"] for r in rows])
    p = np.array([(r["sx"] + r["sy"] + r["sz"]) / 3 for r in rows])
    phig = np.array([r["phi_gp_mean"] for r in rows])
    fsg = np.array([r["phiscale_gp_mean"] for r in rows])
    s = pe > 0.5
    # Gauss point next to each node: local coordinate +-1/sqrt(3) on the node's side
    gp = [{k.strip(): float(v) for k, v in r.items()} for r in csv.DictReader(open(wd / "frozen_gp.csv"))]
    byel = {int(r["elem"]): c for r, c in zip(rows, cen)}
    errs = []
    for r in gp:
        x = np.array([r["x"], r["y"], r["z"]])
        xc = byel[int(r["elem"])]
        g = xc + np.sign(x - xc) * (h / 2) / np.sqrt(3.0)
        # trilinear interpolation of the nodal phi at g
        w = [(1 + np.sign(cc) * (g - xc) / (h / 2)) / 2 for cc in corners]
        ref = sum(np.prod(wi) * seed_phi((xc + cc)[None, :])[0] for wi, cc in zip(w, corners))
        errs.append(abs(r["phi_gp"] - ref))
    out = dict(n=n, E=E, nu=NU, f=FV, k_alpha=KA, T=T, steps=NSTEP,
               kinematics="finite strain (NLGEOM), SOLID185 B-bar; thesis USERMAT core, phi = nodal temperature",
               seed_rule="elements with mean nodal phi > 0.5", seed_elements=int(s.sum()),
               seed_alpha_minus_1=float(xi[s].mean()), seed_vM=float(vm[s].mean()), seed_vM_over_E=float(vm[s].mean() / E),
               max_vM=float(vm.max()), max_vM_centroid=cen[vm.argmax()].round(6).tolist(),
               check_gp_phi_vs_trilinear_max_abs=float(max(errs)), check_gp_points=len(errs),
               check_gp_phi_mean_vs_nodal_mean_max_abs=float(np.abs(phig - pe).max()),
               elements=dict(cx=cen[:, 0].round(6).tolist(), cy=cen[:, 1].round(6).tolist(), cz=cen[:, 2].round(6).tolist(),
                             phi_mean=pe.round(7).tolist(), alpha_minus_1=[float(f"{v:.7g}") for v in xi],
                             vM=[float(f"{v:.7g}") for v in vm], p=[float(f"{v:.7g}") for v in p],
                             phiscale_gp_mean=[float(f"{v:.7g}") for v in fsg]))
    dst = REPO / "ansys_usermat" / "apdl" / "results" / "2026-10-uel" / f"seed_frozen_ansys_{n}.json"
    dst.write_text(json.dumps(out))
    print({k: v for k, v in out.items() if k != "elements"})
    print("wrote", dst)


if __name__ == "__main__":
    main(*sys.argv[1:3])
