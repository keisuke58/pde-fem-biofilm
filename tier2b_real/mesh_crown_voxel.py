"""Robust real-tooth crown volume mesh by VOXELISATION (dependency-free; medical STLs self-intersect
and break tetgen's PLC, so we sample the interior instead).  Reuses the crown extraction/registration
from mesh_crown_real.py, voxel-fills the interior via a +z even-odd ray test against the crown surface,
and splits each interior voxel into 6 C3D4 tets.  Writes cache_crown_real.npz (nodes, tets) so
build_assembly.py job '*crownreal*' uses it.

Blocky surface (validation/Panel-B body, not the smooth Panel-A overlay).  Run in gmsh_env.
  python mesh_crown_voxel.py [voxel_mm=0.35]
"""
import sys
import numpy as np

sys.path.insert(0, "/home/nishioka/IKM_Hiwi/FEM/tier2b_real")
from mesh_crown_real import crown_tris, OUT, X0, Y0, Z_BASE, Z_OCC  # noqa: E402

H = float(sys.argv[1]) if len(sys.argv) > 1 else 0.35


def inside_mask(tris, gx, gy, gz):
    """Even-odd: a voxel centre is inside if a +z ray crosses the surface an odd number of times above
    it.  Per (x,y) column: find triangles whose xy-projection contains (x,y), get the ray-z at the hit,
    then a cell is inside iff an odd number of hits lie above it."""
    v0, v1, v2 = tris[:, 0], tris[:, 1], tris[:, 2]
    x0, y0 = v0[:, 0], v0[:, 1]
    e1x, e1y = v1[:, 0] - x0, v1[:, 1] - y0
    e2x, e2y = v2[:, 0] - x0, v2[:, 1] - y0
    det = e1x * e2y - e2x * e1y
    det = np.where(np.abs(det) < 1e-12, 1e-12, det)
    mask = np.zeros((len(gx), len(gy), len(gz)), bool)
    for i, x in enumerate(gx):
        px = x - x0
        for j, y in enumerate(gy):
            py = y - y0
            u = (px * e2y - e2x * py) / det
            v = (e1x * py - px * e1y) / det
            hit = (u >= 0) & (v >= 0) & (u + v <= 1.0)
            if not hit.any():
                continue
            zc = v0[hit, 2] + u[hit] * (v1[hit, 2] - v0[hit, 2]) + v[hit] * (v2[hit, 2] - v0[hit, 2])
            zc = np.sort(zc)
            # cell inside iff #hits strictly above its centre is odd
            above = len(zc) - np.searchsorted(zc, gz, side="right")
            mask[i, j] = (above % 2) == 1
    return mask


# the 6-tet split of a unit voxel (corner indices 0..7; 0=(0,0,0),1=(1,0,0),2=(1,1,0),3=(0,1,0),
# 4..7 = same with z+1) -- a standard non-degenerate Kuhn decomposition.
VOX = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
CORNER = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                   [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])


def main():
    tris = crown_tris()
    lo = tris.reshape(-1, 3).min(0); hi = tris.reshape(-1, 3).max(0)
    gx = np.arange(lo[0] + H / 2, hi[0], H)
    gy = np.arange(lo[1] + H / 2, hi[1], H)
    gz = np.arange(max(lo[2], Z_BASE) + H / 2, hi[2], H)
    mask = inside_mask(tris, gx, gy, gz)
    print("voxels inside=%d / grid=%dx%dx%d (h=%.2f)" % (mask.sum(), len(gx), len(gy), len(gz), H))
    # build shared corner nodes
    nid = {}; nodes = []
    def cn(ix, iy, iz):
        k = (ix, iy, iz)
        if k not in nid:
            nid[k] = len(nodes)
            nodes.append([lo[0] + ix * H, lo[1] + iy * H, max(lo[2], Z_BASE) + iz * H])
        return nid[k]
    tets = []
    ii, jj, kk = np.where(mask)
    for i, j, k in zip(ii, jj, kk):
        c = [cn(i + dx, j + dy, k + dz) for (dx, dy, dz) in CORNER]
        for a, b, d, e in VOX:
            tets.append([c[a], c[b], c[d], c[e]])
    nodes = np.array(nodes, float); tets = np.array(tets, np.int64)
    # fix winding (positive volume) for Abaqus
    p = nodes[tets]
    v6 = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
    tets[v6 < 0] = tets[v6 < 0][:, [0, 1, 3, 2]]
    print("crown_real(voxel): %d nodes %d tets ; z[%.1f,%.1f] r_max=%.2f"
          % (len(nodes), len(tets), nodes[:, 2].min(), nodes[:, 2].max(),
             np.hypot(nodes[:, 0] - X0, nodes[:, 1] - Y0).max()))
    np.savez(f"{OUT}/cache_crown_real.npz", nodes=nodes, tets=tets)
    print("wrote cache_crown_real.npz (voxel real-tooth crown) -- originals untouched")


if __name__ == "__main__":
    main()
