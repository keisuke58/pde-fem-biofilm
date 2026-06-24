"""Tier-2(b) step 1: volume-mesh the REAL mandible and crop a tet region around teeth 23/24.

The whole Open-Full-Jaw P1 mandible STL (watertight) is turned into a volume mesh, then we keep
only the tetrahedra whose centroid falls inside a crop box around the tooth-23 (implant site) and
tooth-24 (adjacent natural tooth) alveolus. Keeping WHOLE tets gives a valid sub-mesh with a real
(ragged) outer bone surface and the real alveolar sockets preserved. We then test whether the
tooth-24 root region is a VOID (socket) in the bone -- the key assumption for the TIE assembly.

Run (in gmsh_env, with LD_LIBRARY_PATH=$CONDA_PREFIX/lib):
    python mesh_bone_region.py
Writes: bone_region.npz  (nodes, tets[, ] re-indexed), and prints socket diagnostics.
"""
import gmsh
import numpy as np

BASE = "/home/nishioka/IKM_Hiwi/FEM/external_tooth_models/OpenJaw_Dataset/Patient_1"
MANDIBLE = f"{BASE}/Mandible/P1_Mandible.stl"
OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real/bone_region.npz"

# crop box around teeth 23/24 (from STL bbox inspection), with margin
X0, X1 = -76.0, -57.0
Y0, Y1 = -50.0, -34.0
Z0, Z1 = 14.0, 41.0
LC = 1.1          # target tet edge length (mm) -- coarse, refined enough for bone

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.merge(MANDIBLE)

# build a discrete-geometry volume from the STL surface
gmsh.model.mesh.classifySurfaces(40 * np.pi / 180.0, True, True, 40 * np.pi / 180.0)
gmsh.model.mesh.createGeometry()
surfs = gmsh.model.getEntities(2)
loop = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfs])
gmsh.model.geo.addVolume([loop])
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMin", LC)
gmsh.option.setNumber("Mesh.MeshSizeMax", LC)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay
gmsh.model.mesh.generate(3)

# pull tets
etypes, etags, enodes = gmsh.model.mesh.getElements(3)
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
ncoords = np.array(ncoords).reshape(-1, 3)
nid2idx = {int(t): i for i, t in enumerate(ntags)}
tet_conn = None
for et, en in zip(etypes, enodes):
    if et == 4:  # 4-node tetra
        tet_conn = np.array(en, dtype=np.int64).reshape(-1, 4)
assert tet_conn is not None, "no tets generated"
print(f"whole mandible: {len(ncoords)} nodes, {len(tet_conn)} tets")

# centroids
conn_idx = np.vectorize(nid2idx.get)(tet_conn)
cent = ncoords[conn_idx].mean(axis=1)

inbox = ((cent[:, 0] >= X0) & (cent[:, 0] <= X1) &
         (cent[:, 1] >= Y0) & (cent[:, 1] <= Y1) &
         (cent[:, 2] >= Z0) & (cent[:, 2] <= Z1))
sub = conn_idx[inbox]
print(f"cropped to box: {inbox.sum()} tets")

# re-index nodes used by the sub-mesh
used = np.unique(sub)
remap = -np.ones(len(ncoords), dtype=np.int64)
remap[used] = np.arange(len(used))
sub_nodes = ncoords[used]
sub_tets = remap[sub]

# socket diagnostic: are there bone tets in the tooth-24 root region?
# tooth-24 root centroid ~ (-63.9,-41.2,22) (lower z = root)
for label, c in (("tooth24 root", (-63.9, -41.2, 22.0)),
                 ("tooth23 root", (-69.4, -41.0, 22.0)),
                 ("tooth24 crown", (-63.9, -41.2, 36.0))):
    sc = np.array(c)
    d = np.linalg.norm(sub_nodes[sub_tets].mean(axis=1) - sc, axis=1)
    near = (d < 2.0).sum()
    print(f"  bone tets within 2mm of {label} {c}: {near}")

np.savez(OUT, nodes=sub_nodes, tets=sub_tets)
print("wrote", OUT, "nodes=", len(sub_nodes), "tets=", len(sub_tets))
gmsh.finalize()
