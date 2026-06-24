"""Ceramic occlusal crown (body-of-revolution) seated on the generic implant abutment, meshed in the
REAL mandible frame on the tooth-23 implant axis.  Adds the LOAD-BEARING restoration that the bare
screw+abutment model lacked: the occlusal load is no longer applied at the abutment top (z=32.5, moment
arm = 0 above the platform) but transmitted through this crown to its occlusal table at z~40.8, so the
30deg-oblique bite now acts through a realistic crown-height moment arm -> bending of the implant neck
and the peri-implant crestal bone.

A monolithic lithium-disilicate-class ceramic crown (E~95 GPa): cervical emergence flaring from the
1.8 mm abutment radius out to a ~3.8 mm contour, tapering to a domed occlusal table.  Same .npz format
(nodes, tets) as mesh_generic_implant.py; tied to the abutment top in build_assembly.py (job *crown*).

Run in gmsh_env with LD_LIBRARY_PATH=$CONDA_PREFIX/lib:  python mesh_crown.py
"""
import os

import gmsh
import numpy as np

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
X0, Y0 = -69.4, -41.0          # implant / crown axis, real mandible frame (matches mesh_generic_implant)
Z_SEAT = 32.5                  # crown intaglio seat roof = abutment top (Z_ABUT); load-transfer plane
Z_MARGIN = 31.0                # crown cervical margin, ~2 mm above bone crest (29) -> 2 mm metal collar
Z_OCC = 38.0                   # occlusal table = adjacent natural-tooth plane (DENTIN top=38.0)
RAB = 1.8                      # abutment-post radius; crown bore = RAB+gap sheathes the post
RBORE = 1.85                   # crown internal bore radius (0.05 mm cement gap over the post)
LC = float(os.environ.get("CROWN_LC", "0.40"))   # crown tet edge (mm); env-overridable for h-convergence sweeps


def meridian():
    """(r,z) cross-section of a HOLLOW crown that SHEATHES the abutment post (a cap with a central
    socket), one side of the axis.  Outer face: occlusal table -> contour -> cervical margin at z=31
    (so only a ~2 mm transmucosal collar of the Ti abutment shows); inner face: a RBORE socket that
    slips over the abutment post and roofs onto the abutment top at Z_SEAT (the load-transfer seat).
    Clinical crown height ~7 mm, realistic crown:collar proportion."""
    p = [(0.0, Z_OCC),                  # occlusal table centre (domed)
         (2.4, Z_OCC - 0.4),            # occlusal table edge (functional cusp)
         (3.8, Z_SEAT + 2.0),           # height of contour (max bulge)
         (3.0, Z_MARGIN + 0.6),         # axial wall taper toward margin
         (2.5, Z_MARGIN),               # cervical margin (outer, skirt bottom)
         (RBORE, Z_MARGIN),             # margin inner lip (bottom of the socket)
         (RBORE, Z_SEAT),               # bore wall up to the seat plane (sheathes the post)
         (0.0, Z_SEAT)]                 # seat roof on axis (sits on the abutment top)
    return p


def main():
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("crown")
    mer = meridian()
    tags = [gmsh.model.occ.addPoint(X0 + r, Y0, z) for (r, z) in mer]
    lines = [gmsh.model.occ.addLine(tags[i], tags[i + 1]) for i in range(len(tags) - 1)]
    lines.append(gmsh.model.occ.addLine(tags[-1], tags[0]))   # close along the axis
    cl = gmsh.model.occ.addCurveLoop(lines)
    surf = gmsh.model.occ.addPlaneSurface([cl])
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
    print("crown: %d nodes %d tets ; z[%.1f,%.1f] r_max=%.2f"
          % (len(nodes), len(tets), nodes[:, 2].min(), nodes[:, 2].max(),
             np.hypot(nodes[:, 0] - X0, nodes[:, 1] - Y0).max()))
    np.savez(f"{OUT}/cache_crown.npz", nodes=nodes, tets=tets)
    print("wrote cache_crown.npz (ceramic crown) -- originals untouched")


if __name__ == "__main__":
    main()
