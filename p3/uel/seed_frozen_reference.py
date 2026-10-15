"""seed_frozen_reference.py -- the UEL side of a cross-code check with ANSYS (request A-16): the hemispherical
seed of make_seed_inp.py with the phase frozen, so that only Eq. 36 of Klempt et al. 2024 and the
finite-strain mechanics act. Same problem in ANSYS (SOLID185 + the thesis USERMAT, alpha_dot = k_alpha phi):
the two codes share no source.

Problem: unit cube, n^3 hexahedra; nodal phi = seed_phi (0.9 inside radius 0.25 about (0.5, 0.5, 0), tanh edge
of width 0.03, 1e-4 outside), held fixed; growth alpha_dot = k_alpha phi at each Gauss point (phi interpolated
from the nodes), alpha(0) = 1; stiffness (phi^2 + f) E, neo-Hooke deviator + K/2 (J_e - 1)^2, F_e = F/alpha;
u = 0 on z = 0, all other faces traction free; E = 10, nu = 0.49, f = 1e-3, k_alpha = 1e-3, T = 1.1 in 11 steps;
mean dilatation (the analogue of SOLID185's B-bar).

    BIOFILM_ECO_CASE=<dir>/uel_eco.txt python seed_frozen_reference.py [n]     (make_seed_inp.py writes it)
writes results/seed_frozen_<n>.json: per element centroid, alpha - 1 and von Mises (Gauss-point means), the
seed means (elements with mean nodal phi > 0.5) and the largest von Mises.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from driver import NDN, NGP, NSV, props_base  # noqa: E402
from emulate_abaqus import NSVT, F, increment  # noqa: E402
from make_seed_inp import cube, seed_phi  # noqa: E402

E, NU, FV, KA, T, NSTEP = 10.0, 0.49, 1e-3, 1e-3, 1.1, 11


def run(n):
    nodes, elems = cube(n)
    liv = 0.999 * 0.999
    P = props_base(E=E, nu=NU, f=FV, lam=0.0, beta=0.0, L0=KA / (liv * np.log(10.0)), Y=1.0, mu_b=np.log(0.1),
                   s=1e-14, phi_min=2.0)
    props = F(P)
    SV = []
    for _ in elems:                           # one species, living fraction ~ 1, frozen
        sv = np.zeros((NGP, NSV)); sv[:, 0] = 0.999; sv[:, 5] = 0.001; sv[:, 6] = 0.999; sv[:, 13] = 1.0
        SV.append(F(np.r_[sv.ravel(), np.zeros(NSVT - NSV * NGP)]))
    U = np.zeros(NDN * len(nodes)); U[3::NDN] = seed_phi(nodes); U[4::NDN] = 1.0
    fixed = np.zeros(len(U), bool); fixed[3::NDN] = fixed[4::NDN] = True
    base = np.isclose(nodes[:, 2], 0.0)
    for d in range(3):
        fixed[NDN * np.where(base)[0] + d] = True
    for k in range(NSTEP):
        U, SV, it = increment(U.copy(), U, fixed, nodes, elems, SV, props, T / NSTEP, k + 1, ttot=k * T / NSTEP,
                              finite=True, md=True)
    xi = np.array([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv in SV])
    st = np.array([sv[NSV * NGP:].reshape(NGP, 2).mean(0) for sv in SV])
    pe = np.array([seed_phi(nodes[el]).mean() for el in elems])
    cen = np.array([nodes[el].mean(0) for el in elems])
    s = pe > 0.5
    return dict(n=n, E=E, nu=NU, f=FV, k_alpha=KA, T=T, steps=NSTEP, kinematics="finite strain, mean dilatation",
                seed_rule="elements with mean nodal phi > 0.5", seed_elements=int(s.sum()),
                seed_alpha_minus_1=float(xi[s].mean()), seed_vM=float(st[s, 0].mean()), seed_vM_over_E=float(st[s, 0].mean() / E),
                max_vM=float(st[:, 0].max()), max_vM_centroid=cen[st[:, 0].argmax()].tolist(),
                elements=dict(cx=cen[:, 0].round(6).tolist(), cy=cen[:, 1].round(6).tolist(), cz=cen[:, 2].round(6).tolist(),
                              phi_mean=pe.round(7).tolist(), alpha_minus_1=[float(f"{v:.7g}") for v in xi],
                              vM=[float(f"{v:.7g}") for v in st[:, 0]], p=[float(f"{v:.7g}") for v in st[:, 1]]))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    r = run(n)
    out = HERE / "results"; out.mkdir(exist_ok=True)
    (out / f"seed_frozen_{n}.json").write_text(json.dumps(r))
    print({k: v for k, v in r.items() if k != "elements"})
