"""ISO 14801-style implant test coupon for the generic-implant DESIGN studies (thread-dimension sweep,
C3D4-vs-C3D10 accuracy, axial-vs-oblique load).  A parametric titanium screw is osseointegrated
(conforming, bonded by shared mesh) in a bone holder cylinder, embedded to 3 mm below the platform;
the abutment is loaded by ~100 N at a chosen angle, the holder is clamped (ISO 14801 boundary).

This isolates the implant-local physics that the full-mandible tier2b_generic model cannot resolve
cheaply.  gmsh node order for the 10-node tet (type 11) equals the Abaqus C3D10 order, so connectivity
is written directly.

Run (gmsh_env, LD_LIBRARY_PATH=$CONDA_PREFIX/lib):
   python implant_coupon.py  D L pitch taper order angle tag
     D=diam mm  L=body mm  pitch mm  taper(0=cyl,0.3=tapered)  order(1=C3D4,2=C3D10)  angle deg  tag
Writes coupon_<tag>.inp  (+ prints mesh stats).
"""
import sys
import numpy as np
import gmsh

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
D = float(sys.argv[1]) if len(sys.argv) > 1 else 4.1
L = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
PITCH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
TAPER = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
ORDER = int(sys.argv[5]) if len(sys.argv) > 5 else 2
ANGLE = float(sys.argv[6]) if len(sys.argv) > 6 else 30.0
TAG = sys.argv[7] if len(sys.argv) > 7 else "base"

THREAD_DEPTH = 0.40
RB = 4.5                 # bone holder radius
EMBED_EXPOSE = 3.0       # mm of implant coronal exposed (bone top 3 mm below platform), ISO 14801
HAB = 8.0                # abutment height above platform (load point ~11 mm above bone level)
FORCE = 100.0            # N occlusal resultant
ETET, EMAT = ({1: "C3D4", 2: "C3D10"}[ORDER]), None
ELTYPE_GMSH = 4 if ORDER == 1 else 11
MATS = {"TI": (110000., 0.34), "BONE": (13700., 0.30)}


def meridian():
    """(r,z) screw profile: apex(z=0) .. threaded body(z=L) .. abutment top(z=L+HAB)."""
    Rc = D / 2.0 - THREAD_DEPTH
    Ro = D / 2.0
    def taper_f(z):       # core/crest scale: narrower apically when TAPER>0
        t = max(0.0, min(1.0, z / L))
        return 1.0 - TAPER * (1.0 - t)
    p = [(0.0, 0.0), (0.55 * Rc * taper_f(0.6), 0.5)]
    z = 1.0
    p.append((Rc * taper_f(z), z))
    while z + PITCH <= L - 0.3:
        p.append((Ro * taper_f(z + 0.5 * PITCH), z + 0.5 * PITCH))
        p.append((Rc * taper_f(z + PITCH), z + PITCH))
        z += PITCH
    p.append((Ro, L))
    p.append((D / 2.0 * 0.85, L))           # platform shoulder -> abutment
    p.append((D / 2.0 * 0.85, L + HAB))     # abutment wall
    p.append((0.0, L + HAB))                # abutment top on axis
    return p


def main():
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("coupon")
    mer = meridian()
    pts = [gmsh.model.occ.addPoint(r, 0, z) for (r, z) in mer]
    lines = [gmsh.model.occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    lines.append(gmsh.model.occ.addLine(pts[-1], pts[0]))
    cl = gmsh.model.occ.addCurveLoop(lines)
    surf = gmsh.model.occ.addPlaneSurface([cl])
    screw = gmsh.model.occ.revolve([(2, surf)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
    screw_vol = [e for e in screw if e[0] == 3]
    z_bone_top = L - EMBED_EXPOSE
    bone = gmsh.model.occ.addCylinder(0, 0, -1.0, 0, 0, z_bone_top + 1.0, RB)
    # fragment -> conforming shared interface; track which outputs are implant vs bone
    out, omap = gmsh.model.occ.fragment([(3, bone)], screw_vol)
    gmsh.model.occ.synchronize()
    bone_frags = set(t for (d, t) in omap[0] if d == 3)
    screw_frags = set()
    for grp in omap[1:]:
        screw_frags |= set(t for (d, t) in grp if d == 3)
    implant_vols = sorted(screw_frags)
    bone_vols = sorted(bone_frags - screw_frags)
    print("implant vols", implant_vols, "bone vols", bone_vols)

    # mesh-size field: fine on implant, coarser in bone
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.30)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.85)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)   # straight-edge C3D10 (no curved-edge neg-Jacobian)
    gmsh.model.mesh.generate(3)
    if ORDER == 2:
        gmsh.model.mesh.setOrder(2)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid = {int(t): i for i, t in enumerate(ntags)}
    npn = 4 if ORDER == 1 else 10

    def vol_elems(vols):
        out_e = []
        for v in vols:
            et, etags, en = gmsh.model.mesh.getElements(3, v)
            for e, tg, nn in zip(et, etags, en):
                if e == ELTYPE_GMSH:
                    out_e.append(np.array(nn, dtype=np.int64).reshape(-1, npn))
        return np.vstack(out_e) if out_e else np.zeros((0, npn), np.int64)

    imp_e = vol_elems(implant_vols)
    bon_e = vol_elems(bone_vols)
    gmsh.finalize()

    def drop_slivers(conn, vmin=1e-4):
        if len(conn) == 0:
            return conn
        p = np.array([[ncoords[nid[int(c[k])]] for k in range(4)] for c in conn])
        v6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
        keep = np.abs(v6) / 6.0 > vmin
        if (~keep).sum():
            print("  dropped %d sliver elems" % int((~keep).sum()))
        return conn[keep]
    imp_e = drop_slivers(imp_e)
    bon_e = drop_slivers(bon_e)
    print(f"coupon {TAG}: D{D} L{L} p{PITCH} taper{TAPER} {ETET} angle{ANGLE} -> "
          f"implant {len(imp_e)} elems, bone {len(bon_e)} elems, {len(ncoords)} nodes")

    # node sets (coords): clamp = bone bottom + outer cylinder ; load = abutment top on axis.
    # restrict to nodes actually used by 3D elements (revolve leaves orphan geometry nodes).
    used = np.zeros(len(ncoords), dtype=bool)
    for c in np.vstack([imp_e, bon_e]).ravel():
        used[nid[int(c)]] = True
    x, y, z = ncoords[:, 0], ncoords[:, 1], ncoords[:, 2]
    r = np.hypot(x, y)
    clamp = np.where(used & ((z <= -1.0 + 0.4) | (r >= RB - 0.4)))[0]
    ztop = L + HAB
    load = np.where(used & (z >= ztop - 0.4))[0]

    # write INP
    Lout = []; ap = Lout.append
    ap("*HEADING"); ap(" ISO14801-style implant coupon %s" % TAG)
    ap("*NODE")
    for i, (px, py, pz) in enumerate(ncoords, start=1):
        ap(" %d, %.5f, %.5f, %.5f" % (i, px, py, pz))
    gid = 1
    rows = []
    for name, conn in (("TI", imp_e), ("BONE", bon_e)):
        ids = np.arange(gid, gid + len(conn)); gid += len(conn)
        rows.append((name, ids, conn))
        # gmsh tet10 node order differs from Abaqus C3D10 in the last two mid-edge nodes:
        # gmsh 8=mid(2,3), 9=mid(1,3); Abaqus 9=mid(1,3),10=mid(2,3) -> swap cols 8,9.
        perm = list(range(npn)) if ORDER == 1 else [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
        ap("*ELEMENT, TYPE=%s, ELSET=%s" % (ETET, name))
        for e, c in zip(ids, conn):
            cc = [nid[int(c[p])] + 1 for p in perm]
            ap(" " + ",".join(str(v) for v in [e] + cc))
    for m, (E, nu) in MATS.items():
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m); ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))

    def nset(nm, idx):
        idx = np.unique(idx) + 1
        ap("*NSET, NSET=%s" % nm)
        for k in range(0, len(idx), 16):
            ap(" " + ",".join(str(int(v)) for v in idx[k:k + 16]))
    nset("CLAMP", clamp); nset("LOAD", load)

    fx = FORCE * np.sin(np.radians(ANGLE)); fz = -FORCE * np.cos(np.radians(ANGLE))
    nl = max(1, len(np.unique(load)))
    ap("*BOUNDARY"); ap(" CLAMP, 1, 3")
    ap("*STEP, NLGEOM=NO"); ap(" load %g deg" % ANGLE); ap("*STATIC")
    ap("*CLOAD")
    for n in np.unique(load) + 1:
        ap(" %d, 1, %.5f" % (n, fx / nl)); ap(" %d, 3, %.5f" % (n, fz / nl))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    open(f"{OUT}/coupon_{TAG}.inp", "w").write("\n".join(Lout) + "\n")
    import json
    with open(f"{OUT}/coupon_params.jsonl", "a") as f:
        f.write(json.dumps({"tag": TAG, "D": D, "L": L, "pitch": PITCH, "taper": TAPER,
                            "order": ORDER, "angle": ANGLE}) + "\n")
    print("wrote coupon_%s.inp" % TAG)


if __name__ == "__main__":
    main()
