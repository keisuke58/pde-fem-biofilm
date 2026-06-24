"""Peri-implantitis progression coupon: a generic Ti screw in a bone holder with a PARAMETRIC bone
level (marginal bone loss).  As peri-implantitis advances the crestal bone resorbs apically, exposing
more thread and unsupporting the coronal implant.  We lower the bone level in steps and re-solve to
trace the resorption -> stress-redistribution feedback (study (1)); the same meshes feed the
strain-energy mechanostat remodelling (study (2)) and the frictional-contact / micromotion variant
(study (4)).  Biofilm severity (study (3)) enters as the resorption RATE on the time axis, not the mesh.

argv: expose  interface(bond|contact)  order  tag
  expose   = bone level below the platform (mm); larger = more marginal bone loss
  interface= bond (osseointegrated, OCC fragment) | contact (frictional, BIC loss)
  order    = 1 (C3D4) | 2 (C3D10)   [contact uses 1]
Run in gmsh_env (LD_LIBRARY_PATH=$CONDA_PREFIX/lib).  Writes pimp_<tag>.inp.
"""
import sys
import numpy as np
import gmsh

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
EXPOSE = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
IFACE = sys.argv[2] if len(sys.argv) > 2 else "bond"
ORDER = int(sys.argv[3]) if len(sys.argv) > 3 else 2
TAG = sys.argv[4] if len(sys.argv) > 4 else "e3"
# crown lever arm (mm) ABOVE the abutment top: the 30deg bite is applied at a reference point this far
# above the platform-top, kinematically tied to the abutment-top nodes (ISO 14801 offset load). 0 =
# original behaviour (load directly on the abutment top). >0 lengthens the crown-to-implant moment arm.
CROWN_H = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
# eccentric occlusal contact: lateral (+x) offset of the bite point on the crown (mm).  Models a
# functional-cusp / off-axis contact -> an extra bending moment toward +x (study (3), saucerisation).
ECC = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0

D, L, PITCH, THREAD_DEPTH = 4.1, 10.0, 0.8, 0.40
RB, HAB, FORCE, ANGLE = 4.5, 8.0, 100.0, 30.0
MU = 0.3                       # bone-implant friction (contact case)
ETET = {1: "C3D4", 2: "C3D10"}[ORDER]
ELTYPE_GMSH = 4 if ORDER == 1 else 11
MATS = {"TI": (110000., 0.34), "BONE": (13700., 0.30)}
PERM10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]   # gmsh tet10 -> Abaqus C3D10


def meridian():
    Rc, Ro = D / 2.0 - THREAD_DEPTH, D / 2.0
    p = [(0.0, 0.0), (0.55 * Rc, 0.5)]
    z = 1.0
    p.append((Rc, z))
    while z + PITCH <= L - 0.3:
        p.append((Ro, z + 0.5 * PITCH)); p.append((Rc, z + PITCH)); z += PITCH
    p.append((Ro, L)); p.append((D / 2.0 * 0.85, L))
    p.append((D / 2.0 * 0.85, L + HAB)); p.append((0.0, L + HAB))
    return p


def main():
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("pimp")
    mer = meridian()
    pts = [gmsh.model.occ.addPoint(r, 0, z) for (r, z) in mer]
    ln = [gmsh.model.occ.addLine(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    ln.append(gmsh.model.occ.addLine(pts[-1], pts[0]))
    surf = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(ln)])
    screw = [e for e in gmsh.model.occ.revolve([(2, surf)], 0, 0, 0, 0, 0, 1, 2 * np.pi) if e[0] == 3]
    z_bone_top = L - EXPOSE
    bone = gmsh.model.occ.addCylinder(0, 0, -1.0, 0, 0, z_bone_top + 1.0, RB)

    if IFACE == "bond":
        out, omap = gmsh.model.occ.fragment([(3, bone)], screw)
        gmsh.model.occ.synchronize()
        bone_frags = set(t for (d, t) in omap[0] if d == 3)
        screw_frags = set()
        for g in omap[1:]:
            screw_frags |= set(t for (d, t) in g if d == 3)
        imp_vols, bone_vols = sorted(screw_frags), sorted(bone_frags - screw_frags)
    else:  # contact: cut an implant-shaped socket in the bone, keep both as SEPARATE (unmerged)
        # volumes with coincident surfaces -> non-bonded interface for general contact.
        cut, _ = gmsh.model.occ.cut([(3, bone)], screw, removeObject=True, removeTool=False)
        gmsh.model.occ.synchronize()
        imp_vols = [v for (d, v) in screw]
        bone_vols = [v for (d, v) in cut if d == 3]

    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.30)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.85)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.model.mesh.generate(3)
    if ORDER == 2:
        gmsh.model.mesh.setOrder(2)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid = {int(t): i for i, t in enumerate(ntags)}
    npn = 4 if ORDER == 1 else 10

    def vol_elems(vols):
        acc = []
        for v in vols:
            et, _, en = gmsh.model.mesh.getElements(3, v)
            for e, nn in zip(et, en):
                if e == ELTYPE_GMSH:
                    acc.append(np.array(nn, dtype=np.int64).reshape(-1, npn))
        return np.vstack(acc) if acc else np.zeros((0, npn), np.int64)

    imp_e, bon_e = vol_elems(imp_vols), vol_elems(bone_vols)

    # boundary triangles per volume (for contact surfaces); collect now while model is live
    def vol_tris(vols):
        acc = []
        for v in vols:
            for bdim, btag in gmsh.model.getBoundary([(3, v)], oriented=False):
                et, _, en = gmsh.model.mesh.getElements(2, btag)
                for e, nn in zip(et, en):
                    if e in (2, 9):
                        k = 3 if e == 2 else 6
                        acc.append(np.array(nn, dtype=np.int64).reshape(-1, k)[:, :3])
        return np.vstack(acc) if acc else np.zeros((0, 3), np.int64)
    imp_tris = vol_tris(imp_vols) if IFACE == "contact" else None
    bon_tris = vol_tris(bone_vols) if IFACE == "contact" else None
    gmsh.finalize()

    def drop_slivers(conn):
        if len(conn) == 0:
            return conn
        p = np.array([[ncoords[nid[int(c[k])]] for k in range(4)] for c in conn])
        v6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
        return conn[np.abs(v6) / 6.0 > 1e-4]
    imp_e, bon_e = drop_slivers(imp_e), drop_slivers(bon_e)
    print(f"pimp {TAG}: expose={EXPOSE} iface={IFACE} {ETET} -> implant {len(imp_e)} bone {len(bon_e)} "
          f"elems, bone_top z={z_bone_top:.1f}")

    used = np.zeros(len(ncoords), bool)
    for c in np.vstack([imp_e, bon_e]).ravel():
        used[nid[int(c)]] = True
    x, y, z = ncoords[:, 0], ncoords[:, 1], ncoords[:, 2]
    r = np.hypot(x, y)
    clamp = np.where(used & ((z <= -0.6) | (r >= RB - 0.4)))[0]
    load = np.where(used & (z >= L + HAB - 0.4))[0]

    rp = len(ncoords) + 1                       # reference node for the offset crown load (if any)
    Lout = []; ap = Lout.append
    ap("*HEADING"); ap(" peri-implantitis coupon %s expose=%g crown_h=%g" % (TAG, EXPOSE, CROWN_H))
    ap("*NODE")
    for i, (px, py, pz) in enumerate(ncoords, 1):
        ap(" %d, %.5f, %.5f, %.5f" % (i, px, py, pz))
    if CROWN_H > 0:
        ap(" %d, %.5f, 0.0, %.5f" % (rp, ECC, L + HAB + CROWN_H))   # crown occlusal load point (+ECC lateral)
    gid = 1
    for name, conn in (("TI", imp_e), ("BONE", bon_e)):
        ap("*ELEMENT, TYPE=%s, ELSET=%s" % (ETET, name))
        for c in conn:
            cc = [nid[int(c[p])] + 1 for p in (PERM10 if ORDER == 2 else range(npn))]
            ap(" " + ",".join(str(v) for v in [gid] + cc)); gid += 1
    for m, (E, nu) in MATS.items():
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m); ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))

    def nset(nm, idx):
        idx = np.unique(idx) + 1
        ap("*NSET, NSET=%s" % nm)
        for k in range(0, len(idx), 16):
            ap(" " + ",".join(str(int(v)) for v in idx[k:k + 16]))
    nset("CLAMP", clamp); nset("LOAD", load)

    # offset crown load: rigidly tie the abutment-top nodes to the reference point at the crown
    # occlusal height, so the bite is transmitted with the crown-height moment arm (ISO 14801).
    if CROWN_H > 0:
        ap("*SURFACE, TYPE=NODE, NAME=LOADSURF"); ap(" LOAD,")
        ap("*COUPLING, REF NODE=%d, SURFACE=LOADSURF, CONSTRAINT NAME=CROWNLOAD" % rp)
        ap("*KINEMATIC")

    # general contact (Standard): all exterior faces, frictional -- robust, no manual pairing.
    if IFACE == "contact":
        ap("*SURFACE INTERACTION, NAME=BIC"); ap("*FRICTION"); ap(" %.2f" % MU)
        ap("*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=HARD")
        ap("*CONTACT"); ap("*CONTACT INCLUSIONS, ALL EXTERIOR")
        ap("*CONTACT PROPERTY ASSIGNMENT"); ap(" , , BIC")

    fx, fz = FORCE * np.sin(np.radians(ANGLE)), -FORCE * np.cos(np.radians(ANGLE))
    nl = max(1, len(np.unique(load)))
    ap("*BOUNDARY"); ap(" CLAMP, 1, 3")
    step = "*STEP, NLGEOM=NO" if IFACE == "bond" else "*STEP, NLGEOM=YES"
    ap(step); ap(" load")
    if IFACE == "contact":
        ap("*STATIC, STABILIZE=1e-4"); ap(" 0.05, 1.0, 1e-6, 0.2")
    else:
        ap("*STATIC")
    ap("*CLOAD")
    if CROWN_H > 0:                              # whole bite at the offset crown reference point
        ap(" %d, 1, %.5f" % (rp, fx)); ap(" %d, 3, %.5f" % (rp, fz))
    else:
        for n in np.unique(load) + 1:
            ap(" %d, 1, %.5f" % (n, fx / nl)); ap(" %d, 3, %.5f" % (n, fz / nl))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD")
    if IFACE == "contact":
        ap("*CONTACT OUTPUT"); ap(" CSLIP, CPRESS")
    ap("*END STEP")
    open(f"{OUT}/pimp_{TAG}.inp", "w").write("\n".join(Lout) + "\n")
    import json
    open(f"{OUT}/pimp_params.jsonl", "a").write(
        json.dumps({"tag": TAG, "expose": EXPOSE, "iface": IFACE, "order": ORDER,
                    "z_bone_top": z_bone_top, "crown_h": CROWN_H, "ecc": ECC}) + "\n")
    print("wrote pimp_%s.inp" % TAG)


if __name__ == "__main__":
    main()
