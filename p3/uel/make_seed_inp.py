"""make_seed_inp.py -- Abaqus input for the 3D seed of JAXFEM/fuel_front_3d_coupled.py with the unified
user element: unit cube, n^3 elements, a hemispherical seed (radius 0.25, centre (0.5, 0.5, 0)) on the
substrate z = 0, nutrient supplied at the face x = 1.

    python make_seed_inp.py [n] [outdir] [--lam 1 --psiref 0.32] [--antibiotic] [--small] [--nu 0.49] [--md]

Writes uel_seed<n>[_small][_nu49][_md][_ab].inp (the Abaqus job name) and uel_eco.txt (make_column_inp.eco_file, b_i = 1). JPROPS(1) = 1 finite strain
(default), 0 small strain (--small, as the Python 3D model); JPROPS(2) = 1 mean dilatation (--md).
VARIABLES = 128: SVARS(113:128) hold von Mises stress and pressure per Gauss point. With --antibiotic the element
gets PROPS(17) = the total time at which the antibiotic starts (T* = 0.5) and PROPS(18) = alpha* = 0.65. Build and run as in make_column_inp.py (job=uel_seed).
Steps: INIT (1e-8): phi (dof 11) per node as the Python seed, c (dof 12) = c_stop = 0.1, u = 0;
GROW: T = 0.2375 in 60 fixed increments, c = 1 at x = 1, u = 0 at
z = 0. Parameters as fuel_front_3d_coupled.run: D = 0.2, L0 = 2, beta = 0.02, eta_phi = 1, K_phi = 1e-5.
The Python script splits the equations explicitly on a finite-difference grid, so it is a comparison,
not a reference to the last digit (expected within a few per cent at n = 16).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from driver import props_base  # noqa: E402
from make_column_inp import eco_file  # noqa: E402


def cube(n):
    x = np.linspace(0, 1, n + 1)
    nodes = np.array([[a, b, c] for a in x for b in x for c in x])
    nid = lambda i, j, k: (i * (n + 1) + j) * (n + 1) + k
    elems = [[nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
              nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
             for i in range(n) for j in range(n) for k in range(n)]
    return nodes, np.array(elems)


def seed_phi(nodes):
    r = np.sqrt((nodes[:, 0] - 0.5) ** 2 + (nodes[:, 1] - 0.5) ** 2 + nodes[:, 2] ** 2)
    return 0.5 * (1 - np.tanh((r - 0.25) / 0.03)) * 0.9 + 1e-4


def deck(n, path, lam=0.0, psiref=0.0, antibiotic=False, finite=True, nu=0.3, md=False):
    nodes, elems = cube(n)
    P = props_base(s=1.0, D=0.2, L0=2.0, lam=lam, psiref=psiref, nu=nu)
    T, nin, t0 = 0.2375, 60, 1e-8
    if antibiotic:
        P = np.r_[P, t0 + 0.5 * T, 0.65]
    out = ["*HEADING", f"unified two-level model, user element, 3D seed, {n}^3 elements", "*NODE"]
    out += [f"{i + 1}, {x:.10g}, {y:.10g}, {z:.10g}" for i, (x, y, z) in enumerate(nodes)]
    out += [f"*USER ELEMENT, NODES=8, TYPE=U1, PROPERTIES={len(P)}, I PROPERTIES=2, COORDINATES=3, VARIABLES=128", "1, 2, 3, 11, 12",
            "*ELEMENT, TYPE=U1, ELSET=ALLE"]
    out += [f"{e + 1}, " + ", ".join(str(k + 1) for k in el) for e, el in enumerate(elems)]
    vals = [f"{v:.10g}" for v in P]
    vals.append(str(1 if finite else 0))                # JPROPS(1): 1 finite strain, 0 small strain
    vals.append(str(1 if md else 0))                    # JPROPS(2): 1 mean dilatation (nearly incompressible)
    out += ["*UEL PROPERTY, ELSET=ALLE"] + [", ".join(vals[k:k + 8]) for k in range(0, len(vals), 8)]

    def nset(name, mask):
        ids = [str(i + 1) for i in np.where(mask)[0]]
        out.append(f"*NSET, NSET={name}")
        out.extend(", ".join(ids[k:k + 16]) for k in range(0, len(ids), 16))
    nset("ALL", np.ones(len(nodes), bool)); nset("BASE", np.isclose(nodes[:, 2], 0.0))
    nset("SUPPLY", np.isclose(nodes[:, 0], 1.0))
    nset("PROBE", np.isclose(nodes[:, 1], 0.5) & (nodes[:, 2] <= 0.3 + 1e-12))
    phi = seed_phi(nodes)
    out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*STATIC, DIRECT", f"{t0:g}, {t0:g}", "*BOUNDARY", "ALL, 1, 3, 0.0",
            "ALL, 12, 12, 0.1"]
    out += [f"{i + 1}, 11, 11, {p:.10g}" for i, p in enumerate(phi)]
    out += ["*END STEP"]
    out += ["*STEP, NAME=GROW, INC=1000", "*STATIC, DIRECT", f"{T / nin:.10g}, {T:.10g}",
            "*BOUNDARY, OP=NEW", "BASE, 1, 3, 0.0", "SUPPLY, 12, 12, 1.0",
            "*NODE PRINT, NSET=PROBE, FREQUENCY=1000", "U, NT", "*END STEP"]
    path.write_text("\n".join(out) + "\n")
    return nodes, elems, phi


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    get = lambda k, d: float(args[args.index(k) + 1]) if k in args else d
    pos = [a for i, a in enumerate(args) if not a.startswith("--") and (i == 0 or not args[i - 1].startswith("--"))]
    n = int(pos[0]) if pos else 16
    od = Path(pos[1]) if len(pos) > 1 else Path(".")
    od.mkdir(parents=True, exist_ok=True)
    ab = "--antibiotic" in args
    nu, md = get("--nu", 0.3), "--md" in args
    name = (f"uel_seed{n}" + ("_small" if "--small" in args else "") + (f"_nu{int(round(nu * 100))}" if nu != 0.3 else "")
            + ("_md" if md else "") + ("_ab" if ab else "") + ".inp")
    deck(n, od / name, lam=get("--lam", 0.0), psiref=get("--psiref", 0.0), antibiotic=ab,
         finite="--small" not in args, nu=nu, md=md)
    eco_file(od / "uel_eco.txt")
    print(f"wrote {od / name} and {od / 'uel_eco.txt'}")
