"""REAL-tooth-shape ceramic crown as a LOAD-BEARING volume mesh (not just a figure overlay).
Extracts the clinical crown of a real segmented molar STL (OpenJaw P1_Tooth_30), uprights it (long
axis -> z, occlusal up), scales so the occlusal table sits at z=38 (neighbour plane) with the cervical
base at the abutment top z=32.5, caps the open cervical loop -> a watertight surface, then volume-tets
it in gmsh.  Used to VALIDATE the parametric (revolved) crown: peri-implant bone stress is governed by
the load resultant + moment arm, so a real cusped crown should reproduce the ~88 MPa crestal peak.

Writes cache_crown_real.npz (nodes, tets) -- same format as mesh_crown.py; build_assembly.py job
'*crownreal*' uses it.  Run in gmsh_env (LD_LIBRARY_PATH=$CONDA_PREFIX/lib):  python mesh_crown_real.py
"""
import struct
from collections import defaultdict

import numpy as np
import gmsh

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
STL = ("/home/nishioka/IKM_Hiwi/FEM/external_tooth_models/OpenJaw_Dataset/"
       "Patient_1/Teeth/P1_Tooth_30.stl")
X0, Y0 = -69.4, -41.0          # implant axis (real mandible frame)
Z_BASE, Z_OCC = 32.5, 38.0     # cervical base = abutment top ; occlusal table = neighbour plane
CROWN_FRAC = 0.40              # occlusal fraction of the tooth long-axis kept as the clinical crown
LC = 0.40


def load_stl(fn):
    with open(fn, "rb") as f:
        f.read(80); n = struct.unpack("<I", f.read(4))[0]; data = f.read(n * 50)
    rec = np.frombuffer(data, dtype=np.uint8).reshape(n, 50)
    return rec[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(float)


def crown_tris():
    """Real-molar clinical crown, upright (occlusal +z), scaled+registered onto the implant axis,
    cervical loop capped -> watertight triangle soup (n,3,3)."""
    tris = load_stl(STL)
    V = tris.reshape(-1, 3); c = V.mean(0); Vc = V - c
    _, _, vt = np.linalg.svd(Vc - Vc.mean(0), full_matrices=False)
    a = vt[0] / np.linalg.norm(vt[0]); t = Vc @ a
    perp = Vc - np.outer(t, a); rad = np.linalg.norm(perp, axis=1)
    if rad[t > np.percentile(t, 80)].mean() < rad[t < np.percentile(t, 20)].mean():
        a = -a
    zhat = np.array([0, 0, 1.0]); v = np.cross(a, zhat); s = np.linalg.norm(v); cs = a @ zhat
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + vx + vx @ vx * ((1 - cs) / s ** 2)
    P = ((tris.reshape(-1, 3) - c) @ R.T).reshape(-1, 3, 3)
    zc = P[:, :, 2].mean(1); zmax, zmin = zc.max(), zc.min()
    P = P[zc >= zmax - CROWN_FRAC * (zmax - zmin)]                       # occlusal crown band
    sc = (Z_OCC - Z_BASE) / (P[:, :, 2].max() - P[:, :, 2].min()); P = P * sc   # scale to crown height
    m = P.reshape(-1, 3).mean(0)
    P[:, :, 0] += X0 - m[0]; P[:, :, 1] += Y0 - m[1]; P[:, :, 2] += Z_BASE - P[:, :, 2].min()
    P[:, :, 2] = np.where(P[:, :, 2] < Z_BASE + 0.30, Z_BASE, P[:, :, 2])    # snap base rim flat

    # cap the (now planar) cervical loop: walk the boundary loop, triangulate it in 2D (Delaunay,
    # keep only triangles whose centroid lies inside the polygon -> no self-intersection on non-convex).
    from scipy.spatial import Delaunay
    from matplotlib.path import Path as MplPath
    allv = P.reshape(-1, 3); key = np.round(allv, 4)
    uniq, inv = np.unique(key, axis=0, return_inverse=True); tri = inv.reshape(-1, 3)
    ecount = defaultdict(int)
    for f in tri:
        for i in range(3):
            e = tuple(sorted((int(f[i]), int(f[(i + 1) % 3])))); ecount[e] += 1
    bnd = [e for e, n in ecount.items() if n == 1]
    adj = defaultdict(list)
    for u, w in bnd:
        adj[u].append(w); adj[w].append(u)
    start = bnd[0][0]; loop = [start]; prev = -1; cur = start
    while True:
        nxts = [x for x in adj[cur] if x != prev]
        if not nxts:
            break
        prev, cur = cur, nxts[0]
        if cur == start:
            break
        loop.append(cur)
    poly2d = uniq[loop][:, :2]
    dl = Delaunay(poly2d)
    path = MplPath(poly2d)
    cap = []
    for s in dl.simplices:
        cen = poly2d[s].mean(0)
        if path.contains_point(cen):
            p = [[poly2d[s[k], 0], poly2d[s[k], 1], Z_BASE] for k in range(3)]
            cap.append(p)
    print("crown_real: %d side-tris + %d cap-tris ; loop=%d z[%.2f,%.2f]"
          % (len(P), len(cap), len(loop), allv[:, 2].min(), allv[:, 2].max()))
    return np.concatenate([P, np.array(cap)], axis=0)


def write_stl(tris, fn):
    n = len(tris)
    with open(fn, "wb") as f:
        f.write(b"\0" * 80); f.write(struct.pack("<I", n))
        for tr in tris:
            nrm = np.cross(tr[1] - tr[0], tr[2] - tr[0]); nn = np.linalg.norm(nrm)
            nrm = nrm / nn if nn > 0 else nrm
            f.write(struct.pack("<3f", *nrm))
            for vtx in tr:
                f.write(struct.pack("<3f", *vtx))
            f.write(b"\0\0")


def main():
    tris = crown_tris()
    tmp = f"{OUT}/_crown_real_watertight.stl"; write_stl(tris, tmp)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.merge(tmp)
    gmsh.model.mesh.classifySurfaces(40 * np.pi / 180.0, True, True, np.pi)
    gmsh.model.mesh.createGeometry()
    surfs = [e[1] for e in gmsh.model.getEntities(2)]
    loop = gmsh.model.geo.addSurfaceLoop(surfs)
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", LC)
    gmsh.option.setNumber("Mesh.MeshSizeMax", LC)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(3)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    et, _, en = gmsh.model.mesh.getElements(3)
    tet = None
    for e, nn in zip(et, en):
        if e == 4:
            tet = np.vectorize(nid2idx.get)(np.array(nn, dtype=np.int64).reshape(-1, 4))
    gmsh.finalize()
    used = np.unique(tet)
    remap = -np.ones(len(ncoords), dtype=np.int64); remap[used] = np.arange(len(used))
    nodes, tets = ncoords[used], remap[tet]
    print("crown_real volume: %d nodes %d tets ; z[%.1f,%.1f] r_max=%.2f"
          % (len(nodes), len(tets), nodes[:, 2].min(), nodes[:, 2].max(),
             np.hypot(nodes[:, 0] - X0, nodes[:, 1] - Y0).max()))
    np.savez(f"{OUT}/cache_crown_real.npz", nodes=nodes, tets=tets)
    print("wrote cache_crown_real.npz (real-tooth crown volume) -- originals untouched")


if __name__ == "__main__":
    main()
