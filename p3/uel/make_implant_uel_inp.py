"""make_implant_uel_inp.py -- the unified user element on the implant sector of
abaqus_composition/make_implant_inp.py: a 90 degree sector of a biofilm layer around the collar of a
dental implant, hexahedra in (r, theta, z).

    python make_implant_uel_inp.py [outdir] [--nr 2 --nt 6 --nz 8] [--bulge b] [--nu 0.49] [--small] [--no-md]
    BIOFILM_ECO_CASE=<outdir>/uel_eco.txt python make_implant_uel_inp.py [outdir] --run [same options]

Geometry (mm, chosen for this work, as make_implant_inp.py): inner radius 2.05, layer 0.25 thick, 2.0 high.
  - inner face (r = ri): bonded to the titanium, u = 0 (the implant is taken as rigid);
  - theta = 0 and 90 degrees: symmetry planes (u_y = 0 at y = 0, u_x = 0 at x = 0);
  - top face z = 2.0: nutrient c = 1; c = 0.1 elsewhere at the start (c_stop, as make_seed_inp.py);
  - first colonisers: phi = 0.9 on the nodal share of the elements of the inner element layer (as
    make_implant_inp.py), 1e-4 elsewhere; the other faces are free.
  - --bulge b: the inner radius grows to ri + b sin(pi z / 2h) (a crown), as make_implant_inp.py.
Element, properties and the two steps (INIT, GROW: T = 0.2375 in 60 increments) as make_seed_inp.py;
parameters as fuel_front_3d_coupled (D = 0.2, L0 = 2, beta = 0.02).

Without --run it writes uel_implant<nr>x<nt>x<nz>[...].inp, uel_eco.txt and the .json with the seed
elements. With --run it solves the same problem in Python through the Abaqus entry point UEL
(emulate_abaqus.py) and writes results/implant_<nr>x<nt>x<nz>.json: seed means of alpha - 1 and
vM/E (Gauss-point means per element), the largest vM and its centroid.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from driver import props_base  # noqa: E402
from make_column_inp import CTRL, eco_file, overlay  # noqa: E402

RI, THICK, HEIGHT = 2.05, 0.25, 2.0
T, NINC, T0 = 0.2375, 60, 1e-8


def sector(nr, nt, nz, bulge=0.0, ri=RI, thick=THICK, height=HEIGHT):
    """nodes (N, 3), hexahedra (E, 8, 0-based), radial layer index of every element, node (i, j, k) triples"""
    rin = lambda z: ri + bulge * math.sin(0.5 * math.pi * z / height)  # noqa: E731
    nid = lambda i, j, k: i + (nr + 1) * j + (nr + 1) * (nt + 1) * k  # noqa: E731
    nodes, ijk = [], []
    for k in range(nz + 1):
        z = height * k / nz
        for j in range(nt + 1):
            th = 0.5 * math.pi * j / nt
            for i in range(nr + 1):
                r = rin(z) + thick * i / nr
                nodes.append((r * math.cos(th), r * math.sin(th), z))
                ijk.append((i, j, k))
    elems, layer = [], []
    for k in range(nz):
        for j in range(nt):
            for i in range(nr):
                elems.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                              nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
                layer.append(i)
    return np.array(nodes), np.array(elems), np.array(layer), np.array(ijk)


def seed_phi(nodes, elems, layer):
    """nodal share of the inner-layer elements around each node (make_implant_inp.py), 0.9 times it, plus 1e-4"""
    around = np.zeros(len(nodes)); inseed = np.zeros(len(nodes))
    for el, ly in zip(elems, layer):
        around[el] += 1
        if ly == 0:
            inseed[el] += 1
    return 0.9 * inseed / around + 1e-4


def sets(nodes, ijk, tol=1e-9):
    return dict(INNER=ijk[:, 0] == 0, TOP=np.isclose(nodes[:, 2], HEIGHT, atol=tol),
                SYMX=np.abs(nodes[:, 0]) < tol, SYMY=np.abs(nodes[:, 1]) < tol)


def deck(path, nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, name=""):
    P = props_base(s=1.0, D=0.2, L0=2.0, lam=0.0, psiref=0.0, nu=nu)
    S = sets(nodes, ijk)
    out = ["*HEADING", f"unified two-level model, user element, implant sector {name}, {len(elems)} elements", "*NODE"]
    out += [f"{i + 1}, {x:.12g}, {y:.12g}, {z:.12g}" for i, (x, y, z) in enumerate(nodes)]
    out += [f"*USER ELEMENT, NODES=8, TYPE=U1, UNSYMM, PROPERTIES={len(P)}, I PROPERTIES=3, COORDINATES=3, VARIABLES=128",
            "1, 2, 3, 11, 12", "*ELEMENT, TYPE=U1, ELSET=ALLE"]
    out += [f"{e + 1}, " + ", ".join(str(k + 1) for k in el) for e, el in enumerate(elems)]
    vals = [f"{v:.10g}" for v in P] + [str(1 if finite else 0), str(1 if md else 0), "0"]
    out += ["*UEL PROPERTY, ELSET=ALLE"] + [", ".join(vals[k:k + 8]) for k in range(0, len(vals), 8)]
    out += overlay(elems)

    def nset(nm, mask):
        ids = [str(i + 1) for i in np.where(mask)[0]]
        out.append(f"*NSET, NSET={nm}")
        out.extend(", ".join(ids[k:k + 16]) for k in range(0, len(ids), 16))
    nset("ALL", np.ones(len(nodes), bool))
    for k, m in S.items():
        nset(k, m)
    nset("NOTOP", ~S["TOP"])
    phi = seed_phi(nodes, elems, layer)
    out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{T0:g}, {T0:g}", *CTRL, "*BOUNDARY", "ALL, 1, 3, 0.0",
            "NOTOP, 12, 12, 0.1", "TOP, 12, 12, 1.0"]
    out += [f"{i + 1}, 11, 11, {p:.10g}" for i, p in enumerate(phi)]
    out += ["*END STEP"]
    out += ["*STEP, NAME=GROW, INC=1000, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{T / NINC:.10g}, {T:.10g}", *CTRL,
            "*BOUNDARY, OP=NEW", "INNER, 1, 3, 0.0", "SYMX, 1, 1, 0.0", "SYMY, 2, 2, 0.0", "TOP, 12, 12, 1.0", "*END STEP"]
    Path(path).write_text("\n".join(out) + "\n")
    return phi


def run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, nstop=NINC):
    from emulate_abaqus import NDN, NGP, NSV, NSVT, F, increment  # needs BIOFILM_ECO_CASE
    props = F(props_base(s=1.0, D=0.2, L0=2.0, lam=0.0, psiref=0.0, nu=nu))
    S = sets(nodes, ijk)
    SV = [F(np.zeros(NSVT)) for _ in elems]
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy()
    U[3::NDN] = seed_phi(nodes, elems, layer); U[4::NDN] = np.where(S["TOP"], 1.0, 0.1)
    U, SV, _ = increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, T0, 1, finite=finite, md=md)
    fixed = np.zeros(len(U), bool)
    for d in range(3):
        fixed[NDN * np.where(S["INNER"])[0] + d] = True
    fixed[NDN * np.where(S["SYMX"])[0] + 0] = True
    fixed[NDN * np.where(S["SYMY"])[0] + 1] = True
    fixed[NDN * np.where(S["TOP"])[0] + 4] = True
    U[NDN * np.where(S["TOP"])[0] + 4] = 1.0
    for k in range(nstop):
        U, SV, _ = increment(U.copy(), U, fixed, nodes, elems, SV, props, T / NINC, k + 1, ttot=T0 + k * T / NINC,
                             finite=finite, md=md)
    xi = np.array([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv in SV])
    st = np.array([sv[NSV * NGP:].reshape(NGP, 2).mean(0) for sv in SV])
    cen = np.array([nodes[el].mean(0) for el in elems])
    E = 10.0 * 0 + 1.0                                      # props_base default Y = 1: vM/E is vM itself
    s = layer == 0
    im = int(st[:, 0].argmax())
    return dict(seed_elements=int(s.sum()), seed_alpha_minus_1=float(xi[s].mean()), seed_vM_over_E=float(st[s, 0].mean() / E),
                alpha_minus_1_all_mean=float(xi.mean()), max_vM_over_E=float(st[im, 0] / E), max_centroid_xyz=cen[im].tolist(),
                max_centroid_r_theta_z=[float(np.hypot(*cen[im, :2])), float(np.degrees(np.arctan2(cen[im, 1], cen[im, 0]))),
                                        float(cen[im, 2])],
                elements=dict(cx=cen[:, 0].round(6).tolist(), cy=cen[:, 1].round(6).tolist(), cz=cen[:, 2].round(6).tolist(),
                              layer=layer.tolist(), alpha_minus_1=[float(f"{v:.7g}") for v in xi],
                              vM_over_E=[float(f"{v / E:.7g}") for v in st[:, 0]]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", nargs="?", default=".")
    ap.add_argument("--nr", type=int, default=2)
    ap.add_argument("--nt", type=int, default=6)
    ap.add_argument("--nz", type=int, default=8)
    ap.add_argument("--bulge", type=float, default=0.0)
    ap.add_argument("--nu", type=float, default=0.49)
    ap.add_argument("--small", action="store_true")
    ap.add_argument("--no-md", action="store_true")
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    nodes, elems, layer, ijk = sector(a.nr, a.nt, a.nz, a.bulge)
    tag = f"{a.nr}x{a.nt}x{a.nz}" + (f"_b{a.bulge:g}" if a.bulge else "") + ("_small" if a.small else "") + ("_nomd" if a.no_md else "")
    kw = dict(nu=a.nu, md=not a.no_md, finite=not a.small)
    if a.run:
        r = run(nodes, elems, layer, ijk, **kw)
        r.update(nr=a.nr, nt=a.nt, nz=a.nz, bulge=a.bulge, nu=a.nu, kinematics=("small" if a.small else "finite") + " strain"
                 + ("" if a.no_md else ", mean dilatation"), elements_total=len(elems))
        out = HERE / "results"; out.mkdir(exist_ok=True)
        (out / f"implant_{tag}.json").write_text(json.dumps(r))
        print({k: v for k, v in r.items() if k != "elements"})
    else:
        od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
        deck(od / f"uel_implant{tag}.inp", nodes, elems, layer, ijk, name=tag, **kw)
        eco_file(od / "uel_eco.txt")
        (od / f"uel_implant{tag}.json").write_text(json.dumps({"seed": (np.where(layer == 0)[0] + 1).tolist()}))
        print(f"wrote {od / f'uel_implant{tag}.inp'}: {len(elems)} elements, {len(nodes)} nodes")
