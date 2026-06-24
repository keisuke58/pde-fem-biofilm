"""Tier-2(b): FULL real-shape coupled implant + adjacent tooth + alveolar-bone FEM.

A SINGLE conforming tetrahedral mesh of the REAL Open-Full-Jaw P1 mandible segment (cropped around
teeth 23/24) carries embedded multi-materials, each following REAL anatomy:
  * tooth 24 root  -> DENTIN   (point-in-surface test against the real P1_Tooth_24.stl)
  * its ligament   -> PDL      (real P1_PDL_24.stl, minus the dentin it encloses)
  * tooth 23 site  -> TITANIUM IMPLANT (osteotomy cylinder; the natural tooth 23 is "extracted")
  * crestal collar -> BIOFILM  (dysbiotic growth eigenstrain, peri-implant & peri-tooth)
  * remainder      -> BONE     (the solid mandible the STL provides; sockets are CARVED by the
                                reclassification above)
Because it is ONE conforming mesh, load transmits implant <-> shared bone <-> tooth with no TIE/contact
and no gaps -- the genuine Tier-2 coupling on REAL geometry. Two steps: (1) dysbiotic biofilm growth;
(2) occlusal (masticatory) load on both columns.

Units mm, MPa. bone 13.7 GPa, Ti 110 GPa, dentin 18 GPa, PDL 50 MPa, biofilm 1 MPa.
Run (gmsh_env, LD_LIBRARY_PATH=$CONDA_PREFIX/lib):  python build_tier2b.py
Writes tier2b_real.inp + tier2b_meta.npz.
"""
import gmsh
import numpy as np

BASE = "/home/nishioka/IKM_Hiwi/FEM/external_tooth_models/OpenJaw_Dataset/Patient_1"
MANDIBLE = f"{BASE}/Mandible/P1_Mandible.stl"
TOOTH24 = f"{BASE}/Teeth/P1_Tooth_24.stl"
PDL24 = f"{BASE}/PDLs/P1_PDL_24.stl"
OUT_INP = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real/tier2b_real.inp"
OUT_META = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real/tier2b_meta.npz"

# crop box around teeth 23/24
X0, X1 = -76.0, -57.0
Y0, Y1 = -50.0, -34.0
Z0, Z1 = 14.0, 41.0
LC = 1.7                          # bone tet edge (mm)

# implant (osteotomy) at the tooth-23 site, axis along z
IMP_XY = (-69.4, -41.0)
IMP_R = 2.0
IMP_Z = (17.5, 32.5)             # apex .. crest/neck
CREST_Z = 30.0                    # approximate alveolar crest height
EPS_GROWTH = 0.19                 # dysbiotic biofilm growth eigenstrain (DH)

MATS = {"BONE": (13700., 0.30), "TI": (110000., 0.34), "DENTIN": (18000., 0.31),
        "PDL": (50., 0.45), "BIOFILM": (1.0, 0.45)}


def read_tris(path):
    """Return (M,3,3) triangle vertex array from an STL via gmsh."""
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.merge(path)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    etypes, _, enodes = gmsh.model.mesh.getElements(2)
    tri = None
    for et, en in zip(etypes, enodes):
        if et == 2:
            tri = np.array(en, dtype=np.int64).reshape(-1, 3)
    gmsh.finalize()
    idx = np.vectorize(nid2idx.get)(tri)
    return ncoords[idx]


def inside_surface(pts, tris):
    """+x ray-cast parity test: is each point inside the closed triangulated surface?
    pts (N,3), tris (M,3,3).  Returns bool (N,)."""
    res = np.zeros(len(pts), dtype=bool)
    tmin = tris.reshape(-1, 3).min(0); tmax = tris.reshape(-1, 3).max(0)
    inb = ((pts[:, 1] >= tmin[1]) & (pts[:, 1] <= tmax[1]) &
           (pts[:, 2] >= tmin[2]) & (pts[:, 2] <= tmax[2]) & (pts[:, 0] <= tmax[0]))
    cand = np.where(inb)[0]
    if len(cand) == 0:
        return res
    A, B, C = tris[:, 0], tris[:, 1], tris[:, 2]
    Ay, Az = A[:, 1], A[:, 2]; Ax = A[:, 0]
    v0y, v0z = B[:, 1] - Ay, B[:, 2] - Az           # B-A in yz
    v1y, v1z = C[:, 1] - Ay, C[:, 2] - Az           # C-A in yz
    Bx, Cx = B[:, 0], C[:, 0]
    d00 = v0y * v0y + v0z * v0z
    d01 = v0y * v1y + v0z * v1z
    d11 = v1y * v1y + v1z * v1z
    denom = d00 * d11 - d01 * d01
    ok = np.abs(denom) > 1e-14
    Ay, Az, Ax, Bx, Cx = Ay[ok], Az[ok], Ax[ok], Bx[ok], Cx[ok]
    v0y, v0z, v1y, v1z = v0y[ok], v0z[ok], v1y[ok], v1z[ok]
    d00, d01, d11, denom = d00[ok], d01[ok], d11[ok], denom[ok]
    CH = 400
    P = pts[cand]
    counts = np.zeros(len(cand), dtype=np.int64)
    for s in range(0, len(P), CH):
        q = P[s:s + CH]
        py = q[:, 1][:, None] - Ay[None, :]
        pz = q[:, 2][:, None] - Az[None, :]
        d20 = py * v0y[None, :] + pz * v0z[None, :]
        d21 = py * v1y[None, :] + pz * v1z[None, :]
        vv = (d11[None, :] * d20 - d01[None, :] * d21) / denom[None, :]
        ww = (d00[None, :] * d21 - d01[None, :] * d20) / denom[None, :]
        uu = 1.0 - vv - ww
        ins = (uu >= -1e-9) & (vv >= -1e-9) & (ww >= -1e-9)
        xint = uu * Ax[None, :] + vv * Bx[None, :] + ww * Cx[None, :]
        cross = ins & (xint >= q[:, 0][:, None])     # triangle is ahead along +x
        counts[s:s + CH] = cross.sum(axis=1)
    res[cand] = (counts % 2) == 1
    return res


def main():
    # 1) volume-mesh whole mandible, crop tets to box -> master conforming bone mesh
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.merge(MANDIBLE)
    gmsh.model.mesh.classifySurfaces(40 * np.pi / 180, True, True, 40 * np.pi / 180)
    gmsh.model.mesh.createGeometry()
    surfs = gmsh.model.getEntities(2)
    loop = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfs])
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", LC)
    gmsh.option.setNumber("Mesh.MeshSizeMax", LC)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    etypes, _, enodes = gmsh.model.mesh.getElements(3)
    for et, en in zip(etypes, enodes):
        if et == 4:
            tet = np.vectorize(nid2idx.get)(np.array(en, dtype=np.int64).reshape(-1, 4))
    gmsh.finalize()
    cent = ncoords[tet].mean(axis=1)
    inbox = ((cent[:, 0] >= X0) & (cent[:, 0] <= X1) & (cent[:, 1] >= Y0) &
             (cent[:, 1] <= Y1) & (cent[:, 2] >= Z0) & (cent[:, 2] <= Z1))
    tet = tet[inbox]
    used = np.unique(tet)
    remap = -np.ones(len(ncoords), dtype=np.int64); remap[used] = np.arange(len(used))
    nodes = ncoords[used]; tets = remap[tet]
    cent = nodes[tets].mean(axis=1)
    print(f"master bone mesh: {len(nodes)} nodes, {len(tets)} tets")

    # 2) read real inclusion surfaces
    t24 = read_tris(TOOTH24)
    p24 = read_tris(PDL24)

    # 3) classify each tet centroid -> material
    in_t24 = inside_surface(cent, t24)
    in_p24 = inside_surface(cent, p24)
    rxy = np.hypot(cent[:, 0] - IMP_XY[0], cent[:, 1] - IMP_XY[1])
    in_imp = (rxy <= IMP_R) & (cent[:, 2] >= IMP_Z[0]) & (cent[:, 2] <= IMP_Z[1])

    mat = np.array(["BONE"] * len(tets), dtype=object)
    mat[in_p24 & ~in_t24] = "PDL"        # ligament shell (pdl minus dentin)
    mat[in_t24] = "DENTIN"               # tooth root
    mat[in_imp & ~in_t24 & ~in_p24] = "TI"   # implant (tooth-23 site; don't overwrite tooth-24)

    # 4) biofilm collars: thin crestal annulus just outside each column
    rxy_t = np.hypot(cent[:, 0] - (-63.9), cent[:, 1] - (-41.2))   # tooth-24 axis (approx)
    band = (cent[:, 2] >= CREST_Z - 2.5) & (cent[:, 2] <= CREST_Z + 1.0)
    # peri-implant collar: annulus IMP_R..IMP_R+0.9 around implant, currently BONE
    coll_i = band & (rxy > IMP_R) & (rxy <= IMP_R + 0.9) & (mat == "BONE")
    # peri-tooth collar: just outside tooth-24 root (in PDL/bone), near crest
    coll_t = band & (rxy_t > 1.6) & (rxy_t <= 2.6) & np.isin(mat, np.array(["BONE", "PDL"], dtype=object))
    mat[coll_i | coll_t] = "BIOFILM"

    counts = {m: int((mat == m).sum()) for m in MATS}
    print("material tet counts:", counts)

    # 5) node sets: fix the crop-box artificial boundaries (segment far-field); load column tops
    nx, ny, nz = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    tol = 0.6
    fixed = np.where((nz <= Z0 + tol) | (nx <= X0 + tol) | (nx >= X1 - tol) |
                     (ny <= Y0 + tol) | (ny >= Y1 - tol))[0]

    # column tops: highest-z nodes belonging to DENTIN (tooth) and TI (implant) tets
    def top_nodes(material, frac_top=1.2):
        sel = tets[mat == material]
        if len(sel) == 0:
            return np.array([], dtype=np.int64)
        nd = np.unique(sel)
        zmax = nz[nd].max()
        return nd[nz[nd] >= zmax - frac_top]
    too_top = top_nodes("DENTIN")
    imp_top = top_nodes("TI")
    print(f"loaded nodes: tooth-top={len(too_top)} implant-top={len(imp_top)} fixed={len(fixed)}")

    # 6) write Abaqus INP
    L = []; ap = L.append
    ap("*HEADING"); ap(" Tier-2(b) real-shape implant + tooth 24 + mandible, coupled through bone")
    ap("*NODE")
    for i, (x, y, z) in enumerate(nodes, start=1):
        ap(" %d, %.5f, %.5f, %.5f" % (i, x, y, z))
    for m in MATS:
        idx = np.where(mat == m)[0]
        if len(idx) == 0:
            continue
        ap("*ELEMENT, TYPE=C3D4, ELSET=%s" % m)
        for e in idx:
            n = tets[e] + 1
            ap(" %d, %d, %d, %d, %d" % (e + 1, n[0], n[1], n[2], n[3]))
    for m, (E, nu) in MATS.items():
        if (mat == m).sum() == 0:
            continue
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m)
        ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))
        if m == "BIOFILM":
            ap("*EXPANSION"); ap(" %.6f" % EPS_GROWTH)

    def nset(name, ids):
        ids = np.unique(ids) + 1
        ap("*NSET, NSET=%s" % name)
        for k in range(0, len(ids), 16):
            ap(" " + ",".join(str(int(v)) for v in ids[k:k + 16]))
    nset("FIXED", fixed)
    nset("TOOTOP", too_top)
    nset("IMPTOP", imp_top)
    nset("ALLN", np.arange(len(nodes)))

    ap("*INITIAL CONDITIONS, TYPE=TEMPERATURE"); ap(" ALLN, 0.0")
    ap("*BOUNDARY"); ap(" FIXED, 1, 3")
    ap("*STEP"); ap(" 1) dysbiotic biofilm growth"); ap("*STATIC")
    ap("*TEMPERATURE"); ap(" ALLN, 1.0")
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    ap("*STEP"); ap(" 2) occlusal load on both columns -> coupled through shared bone"); ap("*STATIC")
    ap("*CLOAD")
    if len(too_top):
        f = 60.0 / len(too_top)
        for n in np.unique(too_top) + 1:
            ap(" %d, 3, %.4f" % (n, -f)); ap(" %d, 1, %.4f" % (n, 0.2 * f))
    if len(imp_top):
        f = 60.0 / len(imp_top)
        for n in np.unique(imp_top) + 1:
            ap(" %d, 3, %.4f" % (n, -f)); ap(" %d, 1, %.4f" % (n, 0.2 * f))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    open(OUT_INP, "w").write("\n".join(L) + "\n")
    np.savez(OUT_META, nodes=nodes, tets=tets, mat=mat.astype(str), cent=cent)
    print("wrote", OUT_INP)
    print("wrote", OUT_META)


if __name__ == "__main__":
    main()
