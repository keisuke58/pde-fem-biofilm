"""Generic standard screw-form dental implant (parametric), meshed in the REAL mandible frame at the
tooth-23 site, to REPLACE the patient root-analog in the Tier-2(b) coupled model.

A body-of-revolution titanium screw: tapered apex + V-threaded body (concentric thread rings, the
standard axisymmetric-thread idealisation) + transmucosal abutment for occlusal load.  Default
dimensions are a generic Straumann/Nobel-class implant: Ø4.1 mm x 10 mm, 1.0 mm pitch.

Writes cache_implant_generic.npz (nodes, tets) in the SAME format prep_meshes.py produced. The original
root-analog cache_implant.npz and the tier2b_real results are LEFT UNTOUCHED; build_assembly.py is run
with this cache + a new job name (tier2b_generic) to produce the generic-screw variant alongside it.

Run in gmsh_env with LD_LIBRARY_PATH=$CONDA_PREFIX/lib:  python mesh_generic_implant.py
"""
import gmsh
import numpy as np

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
X0, Y0 = -69.4, -41.0          # tooth-23 axis (implant centre), real mandible frame
Z_APEX, Z_PLATFORM = 19.0, 29.0   # implant body spans apex..platform (crest); 10 mm
RC, RO = 1.55, 2.05            # core / thread-crest radius  (Ø4.1 mm)
PITCH = 1.0                    # thread pitch (mm)
RAB, Z_ABUT = 1.8, 32.5        # abutment radius / top (above crest, load application)
LC = 0.45                      # implant tet edge (mm)


def meridian():
    """(r,z) profile points from apex up the threaded flank to the abutment top, on one side of axis."""
    p = [(0.0, Z_APEX), (RC * 0.55, Z_APEX + 0.4)]      # rounded apex tip
    z = Z_APEX + 1.0
    p.append((RC, z))
    while z + PITCH <= Z_PLATFORM - 0.3:
        p.append((RO, z + 0.5 * PITCH))                  # thread crest
        p.append((RC, z + PITCH))                        # thread root
        z += PITCH
    p.append((RO, Z_PLATFORM))                           # platform edge (crest)
    p.append((RAB, Z_PLATFORM))                          # step in to abutment
    p.append((RAB, Z_ABUT))                              # abutment wall
    p.append((0.0, Z_ABUT))                              # abutment top on axis
    return p


def main():
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("implant")
    mer = meridian()
    # build meridian in the x-z plane at y=Y0 (offset r along +x from the axis x=X0)
    tags = [gmsh.model.occ.addPoint(X0 + r, Y0, z) for (r, z) in mer]
    lines = [gmsh.model.occ.addLine(tags[i], tags[i + 1]) for i in range(len(tags) - 1)]
    lines.append(gmsh.model.occ.addLine(tags[-1], tags[0]))   # close along the axis
    cl = gmsh.model.occ.addCurveLoop(lines)
    surf = gmsh.model.occ.addPlaneSurface([cl])
    # revolve full 2*pi about the axis (X0,Y0,*) // z  -> single solid of revolution
    gmsh.model.occ.revolve([(2, surf)], X0, Y0, 0, 0, 0, 1, 2 * np.pi)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", LC)
    gmsh.option.setNumber("Mesh.MeshSizeMax", LC)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = np.array(ncoords).reshape(-1, 3)
    nid2idx = {int(t): i for i, t in enumerate(ntags)}
    et, _, en = gmsh.model.mesh.getElements(3)
    tet = None
    for e, n in zip(et, en):
        if e == 4:
            tet = np.vectorize(nid2idx.get)(np.array(n, dtype=np.int64).reshape(-1, 4))
    gmsh.finalize()
    used = np.unique(tet)
    remap = -np.ones(len(ncoords), dtype=np.int64); remap[used] = np.arange(len(used))
    nodes, tets = ncoords[used], remap[tet]
    print("generic implant: %d nodes %d tets ; z[%.1f,%.1f] r_max=%.2f"
          % (len(nodes), len(tets), nodes[:, 2].min(), nodes[:, 2].max(),
             np.hypot(nodes[:, 0] - X0, nodes[:, 1] - Y0).max()))
    np.savez(f"{OUT}/cache_implant_generic.npz", nodes=nodes, tets=tets)
    print("wrote cache_implant_generic.npz (generic screw) -- originals untouched")


if __name__ == "__main__":
    main()
