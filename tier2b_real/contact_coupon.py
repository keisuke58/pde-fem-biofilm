"""Study (4): loss of osseointegration -> frictional micromotion (Brunski criterion).
The implant and a bone block with an implant-shaped socket are meshed SEPARATELY (coincident but
unmerged surfaces) and assembled into one model with frictional general contact -- i.e. a NON-
osseointegrated implant (0% BIC, the peri-implantitis end state).  Under the ISO 14801 30deg oblique
load the implant micro-moves against the socket; the peak interface slip is compared to the Brunski
~150 um threshold above which fibrous encapsulation (failed re-osseointegration) is expected.

argv: expose mu tag   (expose = bone level below platform, mm; mu = friction)
Run in gmsh_env.  Writes cc_<tag>.inp.
"""
import sys
import numpy as np
import gmsh

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
EXPOSE = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
MU = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
TAG = sys.argv[3] if len(sys.argv) > 3 else "cc"
D, L, PITCH, TD = 4.1, 10.0, 0.8, 0.40
RB, HAB, FORCE, ANGLE = 4.5, 8.0, 100.0, 30.0
Z_BONE_TOP = L - EXPOSE


def meridian():
    Rc, Ro = D / 2 - TD, D / 2
    p = [(0.0, 0.0), (0.55 * Rc, 0.5)]; z = 1.0; p.append((Rc, z))
    while z + PITCH <= L - 0.3:
        p.append((Ro, z + 0.5 * PITCH)); p.append((Rc, z + PITCH)); z += PITCH
    p.append((Ro, L)); p.append((D / 2 * 0.85, L)); p.append((D / 2 * 0.85, L + HAB)); p.append((0.0, L + HAB))
    return p


def screw_dimtags():
    mer = meridian()
    pts = [gmsh.model.occ.addPoint(r, 0, z) for (r, z) in mer]
    ln = [gmsh.model.occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    ln.append(gmsh.model.occ.addLine(pts[-1], pts[0]))
    s = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(ln)])
    return [e for e in gmsh.model.occ.revolve([(2, s)], 0, 0, 0, 0, 0, 1, 2 * np.pi) if e[0] == 3]


def mesh_part(kind):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(kind)
    sc = screw_dimtags()
    if kind == "implant":
        keep = sc
    else:  # bone block with implant-shaped socket
        bone = gmsh.model.occ.addCylinder(0, 0, -1.0, 0, 0, Z_BONE_TOP + 1.0, RB)
        out, _ = gmsh.model.occ.cut([(3, bone)], sc, removeObject=True, removeTool=True)
        keep = [e for e in out if e[0] == 3]
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.35)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.9)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(3)
    nt, nc, _ = gmsh.model.mesh.getNodes(); nc = np.array(nc).reshape(-1, 3)
    nid = {int(t): i for i, t in enumerate(nt)}
    et, _, en = gmsh.model.mesh.getElements(3)
    tet = None
    for e, nn in zip(et, en):
        if e == 4:
            tet = np.vectorize(nid.get)(np.array(nn, dtype=np.int64).reshape(-1, 4))
    gmsh.finalize()
    used = np.unique(tet); remap = -np.ones(len(nc), np.int64); remap[used] = np.arange(len(used))
    return nc[used], remap[tet]


def clean(conn, nodes):
    p = nodes[conn]
    v6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    out = conn.copy(); neg = v6 < 0; out[neg] = out[neg][:, [0, 1, 3, 2]]
    return out[np.abs(v6) / 6 > 1e-4]


def main():
    ni, ti = mesh_part("implant")
    nb, tb = mesh_part("bone")
    ti = clean(ti, ni); tb = clean(tb, nb)
    nodes = np.vstack([ni, nb]); tb_g = tb + len(ni)
    print(f"contact coupon {TAG}: implant {len(ti)} tets, bone {len(tb)} tets, {len(nodes)} nodes")

    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]; r = np.hypot(x, y)
    nI = len(ni)
    clamp = np.where((np.arange(len(nodes)) >= nI) & ((z <= -0.6) | (r >= RB - 0.4)))[0]
    load = np.where((np.arange(len(nodes)) < nI) & (z >= L + HAB - 0.4))[0]

    Lo = []; ap = Lo.append
    ap("*HEADING"); ap(" contact coupon %s (0%% BIC, mu=%g)" % (TAG, MU))
    ap("*NODE")
    for i, (px, py, pz) in enumerate(nodes, 1):
        ap(" %d, %.5f, %.5f, %.5f" % (i, px, py, pz))
    gid = 1
    for name, conn in (("TI", ti), ("BONE", tb_g)):
        ap("*ELEMENT, TYPE=C3D4, ELSET=%s" % name)
        for c in conn:
            ap(" %d, %d, %d, %d, %d" % (gid, c[0] + 1, c[1] + 1, c[2] + 1, c[3] + 1)); gid += 1
    for m, (E, nu) in (("TI", (110000., 0.34)), ("BONE", (13700., 0.30))):
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m); ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))

    def nset(nm, idx):
        idx = np.unique(idx) + 1; ap("*NSET, NSET=%s" % nm)
        for k in range(0, len(idx), 16):
            ap(" " + ",".join(str(int(v)) for v in idx[k:k + 16]))
    nset("CLAMP", clamp); nset("LOAD", load)

    ap("*SURFACE INTERACTION, NAME=BIC"); ap("*FRICTION"); ap(" %.2f" % MU)
    ap("*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=HARD")
    ap("*CONTACT"); ap("*CONTACT INCLUSIONS, ALL EXTERIOR")
    ap("*CONTACT PROPERTY ASSIGNMENT"); ap(" , , BIC")

    fx, fz = FORCE * np.sin(np.radians(ANGLE)), -FORCE * np.cos(np.radians(ANGLE))
    nl = max(1, len(np.unique(load)))
    ap("*BOUNDARY"); ap(" CLAMP, 1, 3")
    ap("*STEP, NLGEOM=YES"); ap(" oblique load on debonded implant")
    ap("*STATIC, STABILIZE=1e-4"); ap(" 0.05, 1.0, 1e-7, 0.2")
    ap("*CLOAD")
    for n in np.unique(load) + 1:
        ap(" %d, 1, %.5f" % (n, fx / nl)); ap(" %d, 3, %.5f" % (n, fz / nl))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    open(f"{OUT}/cc_{TAG}.inp", "w").write("\n".join(Lo) + "\n")
    print("wrote cc_%s.inp" % TAG)


if __name__ == "__main__":
    main()
