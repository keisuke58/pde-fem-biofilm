"""make_implant_uel_inp.py -- the unified user element on the implant sector of
abaqus_composition/make_implant_inp.py: a 90 degree sector of a biofilm layer around the collar of a
dental implant, hexahedra in (r, theta, z).

    python make_implant_uel_inp.py [outdir] [--nr 2 --nt 6 --nz 8] [--bulge b] [--nu 0.49] [--small] [--no-md] [--seed-layers n]
    BIOFILM_ECO_CASE=<outdir>/uel_eco.txt python make_implant_uel_inp.py [outdir] --run [same options]

--seed-layers n seeds the inner n element layers instead of one (the seed of one element layer thins with
every refinement; n = nr/2 keeps it 0.125 mm thick at every level, for the mesh convergence; tag _sl<n>).

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
(emulate_abaqus.py) and writes results/implant_<nr>x<nt>x<nz>[_b1].json: seed means of alpha - 1 and
vM/E (Gauss-point means per element), the same over the seed interior (without the element rows at z = 0,
z = 2.0 and at the two sector cuts, where the fixed inner face meets a free face), the largest vM/E and
its centroid (over all elements and over the interior), the seed mean pressure (positive in compression) and
the wall time. `implant_convergence.py` tabulates the levels; `compare_implant_gp.py` compares an Abaqus run's
uel_gp.txt with the .json.

Homeostatic-pressure law (16 Oct 2026; unified_uel.f90 module header, homeostatic_test.py): `--p-ref x` sets
PROPS(19) = p_ref / E (the deck's E is 1) and `--p-form` its form (0: p_h = p_ref, as keio_wp2 "E"); tag `_ph<x>`.
Boundary conditions beyond make_implant_inp.py: `--outer fixed` (u = 0 on r = ri + 0.25, a rigid confinement; tag
`_of`), `--bottom uz|fixed` (tags `_bz`, `_bf`), for the boundary-condition dependence with and without the law
(keio_wp2/README.ja.md stage 6-7). `--ninc n` changes the increments of GROW (tag `_n<n>`), a time-step check.
`homeostatic_implant.py` runs and tabulates the matrix. `--nut outer` (tag `_no`) supplies the nutrient on the outer
face r = ri + 0.25 instead of the top (saliva on a crown, as the keio_wp2 tooth of stage 7); `cervical_study.py`
uses it for the cervical line of the tooth-like sector (P3 (e)).
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
# *CONTROLS, PARAMETERS=TIME INCREMENTATION: I0, IR as CTRL, I_C = 30 iterations, I_A = 5 cutbacks per increment
CTRL_AUTO = [*CTRL[:2], "*CONTROLS, PARAMETERS=TIME INCREMENTATION", "20, 30, , 30, , , , 5"]


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


def seed_phi(nodes, elems, layer, seed_layers=1):
    """nodal share of the seeded elements (the inner seed_layers element layers; make_implant_inp.py: one) around
    each node, 0.9 times it, plus 1e-4"""
    around = np.zeros(len(nodes)); inseed = np.zeros(len(nodes))
    for el, ly in zip(elems, layer):
        around[el] += 1
        if ly < seed_layers:
            inseed[el] += 1
    return 0.9 * inseed / around + 1e-4


def sets(nodes, ijk, tol=1e-9):
    return dict(INNER=ijk[:, 0] == 0, OUTER=ijk[:, 0] == ijk[:, 0].max(), TOP=np.isclose(nodes[:, 2], HEIGHT, atol=tol),
                BOTTOM=np.isclose(nodes[:, 2], 0.0, atol=tol), SYMX=np.abs(nodes[:, 0]) < tol, SYMY=np.abs(nodes[:, 1]) < tol)


def element_props(nu, p_ref=0.0, p_form=0):
    """the 16 properties, plus (17, 18: antibiotic off; 19, 20) when the homeostatic-pressure law is on"""
    P = props_base(s=1.0, D=0.2, L0=2.0, lam=0.0, psiref=0.0, nu=nu)
    if p_ref > 0:
        P = np.r_[P, 1e30, 0.0, p_ref, float(p_form)]
    return P


def extra_bc(outer, bottom):
    """boundary lines of the GROW step for the outer face (free | fixed: u = 0) and the bottom face (free | uz:
    u_z = 0 | fixed: u = 0); the defaults are the free faces of make_implant_inp.py"""
    out = []
    if outer == "fixed":
        out.append("OUTER, 1, 3, 0.0")
    if bottom == "uz":
        out.append("BOTTOM, 3, 3, 0.0")
    elif bottom == "fixed":
        out.append("BOTTOM, 1, 3, 0.0")
    return out


def deck(path, nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=1, name="", p_ref=0.0, p_form=0,
         outer="free", bottom="free", nut="top"):
    P = element_props(nu, p_ref, p_form)
    S = sets(nodes, ijk)
    NUT = "TOP" if nut == "top" else "OUTER"
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
    nset("NO" + NUT, ~S[NUT])
    phi = seed_phi(nodes, elems, layer, seed_layers)
    out += ["*STEP, NAME=INIT, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT", f"{T0:g}, {T0:g}", *CTRL, "*BOUNDARY", "ALL, 1, 3, 0.0",
            f"NO{NUT}, 12, 12, 0.1", f"{NUT}, 12, 12, 1.0"]
    out += [f"{i + 1}, 11, 11, {p:.10g}" for i, p in enumerate(phi)]
    out += ["*END STEP"]
    # GROW: automatic incrementation as emulate_abaqus.march: DELTMX (never reached, phi and c are in [0, 1]) switches
    # it on, the data line gives dt0 = DTMAX = T/NINC and DTMIN = T/NINC/4^5; the controls allow 30 iterations (I_C)
    # and 5 cutbacks (I_A) per increment; the cutback factor 0.25 and the growth factor 1.5 are Abaqus's defaults
    out += ["*STEP, NAME=GROW, INC=1000, AMPLITUDE=STEP", "*COUPLED TEMPERATURE-DISPLACEMENT, DELTMX=10",
            f"{T / NINC:.10g}, {T:.10g}, {T / NINC / 4 ** 5:.6g}, {T / NINC:.10g}", *CTRL_AUTO,
            "*BOUNDARY, OP=NEW", "INNER, 1, 3, 0.0", "SYMX, 1, 1, 0.0", "SYMY, 2, 2, 0.0", *extra_bc(outer, bottom),
            f"{NUT}, 12, 12, 1.0", "*END STEP"]
    Path(path).write_text("\n".join(out) + "\n")
    return phi


def element_jk(e, nr, nt):
    """(j, k) of element e: sector() numbers the elements i fastest, then j (theta), then k (z)"""
    return (e // nr) % nt, e // (nr * nt)


def run(nodes, elems, layer, ijk, nu=0.49, md=True, finite=True, seed_layers=1, nstop=NINC, maxit=30, p_ref=0.0, p_form=0,
        outer="free", bottom="free", ninc=None, nut="top"):
    """INIT and GROW in the emulation; GROW with automatic incrementation (emulate_abaqus.march, as the deck).
    nstop: GROW ends at nstop * T/NINC (the full step by default); maxit: Newton iterations per attempt (I_C);
    p_ref, p_form: the homeostatic-pressure law (PROPS(19), PROPS(20) of unified_uel.f90; 0 = off);
    outer, bottom: extra_bc; ninc: increments of the full GROW step (default NINC; a time-step check);
    nut: the face with c = 1, "top" (the implant deck) or "outer" (saliva on the crown, keio_wp2 stage 7)."""
    import time as _t
    from emulate_abaqus import NDN, NGP, NSV, NSVT, F, increment, march  # needs BIOFILM_ECO_CASE
    t_start = _t.time()
    NI = NINC if ninc is None else ninc
    props = F(element_props(nu, p_ref, p_form))
    S = sets(nodes, ijk)
    NUT = "TOP" if nut == "top" else "OUTER"
    SV = [F(np.zeros(NSVT)) for _ in elems]
    U0 = np.zeros(NDN * len(nodes)); U = U0.copy()
    U[3::NDN] = seed_phi(nodes, elems, layer, seed_layers); U[4::NDN] = np.where(S[NUT], 1.0, 0.1)
    U, SV, _ = increment(U0, U, np.ones(len(U), bool), nodes, elems, SV, props, T0, 1, finite=finite, md=md)
    fixed = np.zeros(len(U), bool)
    for d in range(3):
        fixed[NDN * np.where(S["INNER"])[0] + d] = True
    fixed[NDN * np.where(S["SYMX"])[0] + 0] = True
    fixed[NDN * np.where(S["SYMY"])[0] + 1] = True
    if outer == "fixed":
        for d in range(3):
            fixed[NDN * np.where(S["OUTER"])[0] + d] = True
    if bottom in ("uz", "fixed"):
        for d in ((2,) if bottom == "uz" else (0, 1, 2)):
            fixed[NDN * np.where(S["BOTTOM"])[0] + d] = True
    fixed[NDN * np.where(S[NUT])[0] + 4] = True
    U[NDN * np.where(S[NUT])[0] + 4] = 1.0
    log = []
    U, SV, info = march(U, fixed, nodes, elems, SV, props, T if nstop == NINC else nstop * (T / NI), nstop if nstop != NINC else NI,
                        t0=T0, finite=finite, md=md, maxit=maxit, log=log)
    xi = np.array([sv[:NSV * NGP].reshape(NGP, NSV)[:, 12].mean() for sv in SV])
    st = np.array([sv[NSV * NGP:].reshape(NGP, 2).mean(0) for sv in SV])
    cen = np.array([nodes[el].mean(0) for el in elems])
    E = 10.0 * 0 + 1.0                                      # props_base default Y = 1: vM/E is vM itself
    s = layer < seed_layers
    nr, nt, nz = int(layer.max()) + 1, int(ijk[:, 1].max()), int(ijk[:, 2].max())
    jk = np.array([element_jk(e, nr, nt) for e in range(len(elems))])
    # interior: not in the top or bottom row of elements (z = 0, z = HEIGHT) and not in the rows at the
    # sector cuts (theta = 0, 90 degrees), where the fixed inner face meets a free face (possible singularities)
    inner = (jk[:, 0] > 0) & (jk[:, 0] < nt - 1) & (jk[:, 1] > 0) & (jk[:, 1] < nz - 1)
    si = s & inner
    rtz = lambda c: [float(np.hypot(*c[:2])), float(np.degrees(np.arctan2(c[1], c[0]))), float(c[2])]  # noqa: E731
    im = int(st[:, 0].argmax())
    ii = int(np.where(inner)[0][st[inner, 0].argmax()]) if inner.any() else im
    return dict(seed_layers=int(seed_layers), seed_thickness=float(seed_layers * THICK / nr), seed_elements=int(s.sum()), seed_alpha_minus_1=float(xi[s].mean()), seed_vM_over_E=float(st[s, 0].mean() / E),
                seed_interior_elements=int(si.sum()),
                seed_interior_alpha_minus_1=float(xi[si].mean()) if si.any() else None,
                seed_interior_vM_over_E=float(st[si, 0].mean() / E) if si.any() else None,
                seed_p_over_E=float(st[s, 1].mean() / E), seed_max_p_over_E=float(st[s, 1].max() / E),
                alpha_minus_1_all_mean=float(xi.mean()), max_vM_over_E=float(st[im, 0] / E), max_centroid_xyz=cen[im].tolist(),
                max_centroid_r_theta_z=rtz(cen[im]),
                max_interior_vM_over_E=float(st[ii, 0] / E), max_interior_centroid_r_theta_z=rtz(cen[ii]),
                p_ref=float(p_ref), p_form=int(p_form), outer=outer, bottom=bottom, nut=nut, increments_nominal=int(NI),
                wall_time_s=float(_t.time() - t_start),
                incrementation=dict(increments=info["increments"], cutbacks=info["cutbacks"], iterations=info["iterations"],
                                    dt_min=info["dt_min"], dt_base=T / NI,
                                    log=[[k, float(f"{t:.10g}"), float(f"{d:.10g}"), it] for k, t, d, it in log]),
                elements=dict(cx=cen[:, 0].round(6).tolist(), cy=cen[:, 1].round(6).tolist(), cz=cen[:, 2].round(6).tolist(),
                              layer=layer.tolist(), interior=inner.astype(int).tolist(), alpha_minus_1=[float(f"{v:.7g}") for v in xi],
                              vM_over_E=[float(f"{v / E:.7g}") for v in st[:, 0]], p_over_E=[float(f"{v / E:.7g}") for v in st[:, 1]]))


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
    ap.add_argument("--seed-layers", type=int, default=1, help="seed = the inner n element layers (n = nr/2: 0.125 mm at every level)")
    ap.add_argument("--p-ref", type=float, default=0.0, help="homeostatic-pressure law: p_ref / E (PROPS(19)); 0 = off")
    ap.add_argument("--p-form", type=int, default=0, help="its form (PROPS(20)): 0 p_h = p_ref, 1 p_h = p_ref (phi^2 + f), 2 hyperbolic")
    ap.add_argument("--outer", choices=("free", "fixed"), default="free", help="outer face r = ri + 0.25: free or u = 0")
    ap.add_argument("--bottom", choices=("free", "uz", "fixed"), default="free", help="bottom face z = 0: free, u_z = 0 or u = 0")
    ap.add_argument("--ninc", type=int, default=None, help="increments of the GROW step (default 60; a time-step check)")
    ap.add_argument("--nut", choices=("top", "outer"), default="top", help="face with c = 1: top (implant collar) or outer (saliva on a crown)")
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    nodes, elems, layer, ijk = sector(a.nr, a.nt, a.nz, a.bulge)
    tag = f"{a.nr}x{a.nt}x{a.nz}" + (f"_b{a.bulge:g}" if a.bulge else "") + ("_small" if a.small else "") + ("_nomd" if a.no_md else "") + (f"_sl{a.seed_layers}" if a.seed_layers != 1 else "")
    tag += (f"_ph{a.p_ref:g}" + (f"f{a.p_form}" if a.p_form else "") if a.p_ref > 0 else "") + ("_of" if a.outer == "fixed" else "")
    tag += {"free": "", "uz": "_bz", "fixed": "_bf"}[a.bottom] + ("_no" if a.nut == "outer" else "") + (f"_n{a.ninc}" if a.ninc else "")
    kw = dict(nu=a.nu, md=not a.no_md, finite=not a.small, seed_layers=a.seed_layers, p_ref=a.p_ref, p_form=a.p_form,
              outer=a.outer, bottom=a.bottom, nut=a.nut)
    if a.run:
        r = run(nodes, elems, layer, ijk, ninc=a.ninc, **kw)
        r.update(nr=a.nr, nt=a.nt, nz=a.nz, bulge=a.bulge, nu=a.nu, kinematics=("small" if a.small else "finite") + " strain"
                 + ("" if a.no_md else ", mean dilatation"), elements_total=len(elems))
        out = HERE / "results"; out.mkdir(exist_ok=True)
        (out / f"implant_{tag}.json").write_text(json.dumps(r))
        print({k: (v if k != "incrementation" else {kk: vv for kk, vv in v.items() if kk != "log"}) for k, v in r.items() if k != "elements"})
    else:
        od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
        deck(od / f"uel_implant{tag}.inp", nodes, elems, layer, ijk, name=tag, **kw)
        eco_file(od / "uel_eco.txt")
        (od / f"uel_implant{tag}.json").write_text(json.dumps({"seed": (np.where(layer < a.seed_layers)[0] + 1).tolist()}))
        print(f"wrote {od / f'uel_implant{tag}.inp'}: {len(elems)} elements, {len(nodes)} nodes")
