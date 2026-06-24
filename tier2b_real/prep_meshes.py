"""Tier-2(b) mesh prep: volume-mesh the real BONE (mandible crop), DENTIN (tooth 24) and IMPLANT
(tooth 23, root-form) parts and cache them, plus the tooth-24 boundary surface for building the
conforming PDL offset layer downstream.  Run in gmsh_env (LD_LIBRARY_PATH=$CONDA_PREFIX/lib).
"""
import gmsh
import numpy as np

BASE = "/home/nishioka/IKM_Hiwi/FEM/external_tooth_models/OpenJaw_Dataset/Patient_1"
OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"

# tight crop around the tooth 23/24 alveolus (keeps bone-tet count tractable; bone crest ~z29)
X0, X1 = -76.0, -59.0
Y0, Y1 = -47.0, -37.5
Z0, Z1 = 15.0, 31.0


def _coarsen_opts(lc, remesh_surface=True):
    # remesh_surface=True (small clean parts): fully control size, remesh surfaces.
    # remesh_surface=False (complex mandible): keep STL surface triangulation; only coarsen interior.
    if remesh_surface:
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0 if not remesh_surface else lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)


def _geom_from_stl(path):
    gmsh.merge(path)
    gmsh.model.mesh.classifySurfaces(40 * np.pi / 180, True, True, 40 * np.pi / 180)
    gmsh.model.mesh.createGeometry()
    s = gmsh.model.getEntities(2)
    loop = gmsh.model.geo.addSurfaceLoop([e[1] for e in s])
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()


def _get_tets():
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    et, _, en = gmsh.model.mesh.getElements(3)
    tet = None
    for e, n in zip(et, en):
        if e == 4:
            tet = np.vectorize(nid2idx.get)(np.array(n, dtype=np.int64).reshape(-1, 4))
    return ncoords, tet


def _get_surface(ncoords, nid2idx):
    """boundary triangles of the current volume mesh, re-indexed to a compact surface vert array."""
    et, _, en = gmsh.model.mesh.getElements(2)
    tri = None
    for e, n in zip(et, en):
        if e == 2:
            tri = np.vectorize(nid2idx.get)(np.array(n, dtype=np.int64).reshape(-1, 3))
    used = np.unique(tri)
    remap = -np.ones(len(ncoords), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return ncoords[used], remap[tri], used


def mesh_solid(path, lc, want_surface=False):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    _geom_from_stl(path)
    _coarsen_opts(lc)
    gmsh.model.mesh.generate(3)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    et, _, en = gmsh.model.mesh.getElements(3)
    for e, n in zip(et, en):
        if e == 4:
            tet = np.vectorize(nid2idx.get)(np.array(n, dtype=np.int64).reshape(-1, 4))
    surf = _get_surface(ncoords, nid2idx) if want_surface else None
    gmsh.finalize()
    # compact to used nodes
    used = np.unique(tet)
    remap = -np.ones(len(ncoords), dtype=np.int64); remap[used] = np.arange(len(used))
    return ncoords[used], remap[tet], surf


def mesh_bone(lc):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    _geom_from_stl(f"{BASE}/Mandible/P1_Mandible.stl")
    # working config: keep the STL surface triangulation (fine); volume follows it.
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)
    ncoords, tet = _get_tets()
    gmsh.finalize()
    cent = ncoords[tet].mean(axis=1)
    inbox = ((cent[:, 0] >= X0) & (cent[:, 0] <= X1) & (cent[:, 1] >= Y0) &
             (cent[:, 1] <= Y1) & (cent[:, 2] >= Z0) & (cent[:, 2] <= Z1))
    tet = tet[inbox]
    used = np.unique(tet)
    remap = -np.ones(len(ncoords), dtype=np.int64); remap[used] = np.arange(len(used))
    return ncoords[used], remap[tet]


if __name__ == "__main__":
    bn, bt = mesh_bone(1.5)
    print(f"BONE: {len(bn)} nodes {len(bt)} tets")
    np.savez(f"{OUT}/cache_bone.npz", nodes=bn, tets=bt)

    dn, dt, dsurf = mesh_solid(f"{BASE}/Teeth/P1_Tooth_24.stl", 1.0, want_surface=True)
    sverts, sfaces, _ = dsurf
    print(f"DENTIN(tooth24): {len(dn)} nodes {len(dt)} tets ; surface {len(sverts)} v {len(sfaces)} f")
    np.savez(f"{OUT}/cache_dentin.npz", nodes=dn, tets=dt, sverts=sverts, sfaces=sfaces)

    inp, it, _ = mesh_solid(f"{BASE}/Teeth/P1_Tooth_23.stl", 1.0)
    print(f"IMPLANT(tooth23): {len(inp)} nodes {len(it)} tets")
    np.savez(f"{OUT}/cache_implant.npz", nodes=inp, tets=it)
    print("cached all parts")
