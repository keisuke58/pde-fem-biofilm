"""fig_matrix_structure.py -- what is solved where in the coupled element (colloquium backup).

(a) Global stiffness matrix: only the displacements (3 per node), 4^3 hexahedra (125 nodes), trilinear coupling.
(b) Field update at the Gauss points: phi_3D, c, a are updated explicitly; each point uses its 30 nearest
    neighbours (NEM stencil), so this is a sparse matrix-vector product, not a solve.
(c) Point model at one Gauss point: the 12 x 12 Newton Jacobian (5 species: phi_1..5, phi_0, psi_1..5, gamma),
    solved alone at every point.

    python ansys_usermat/apdl/fig_matrix_structure.py -> assets/fig_matrix_structure.pdf/.png
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import figstyle  # noqa: E402
figstyle.apply()

ne = 4; nn = ne + 1
nid = lambda i, j, k: i + nn * (j + nn * k)
N = nn ** 3
K = np.zeros((3 * N, 3 * N), bool)
for i in range(ne):
    for j in range(ne):
        for k in range(ne):
            nodes = [nid(i + a, j + b, k + c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
            dofs = [3 * n + d for n in nodes for d in range(3)]
            K[np.ix_(dofs, dofs)] = True

# Gauss points: 2x2x2 per element, 30 nearest neighbours
g1 = (np.arange(2 * ne) + 0.5) / (2 * ne)
X = np.array([(x, y, z) for z in g1 for y in g1 for x in g1])
D = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)
L = np.zeros_like(D, bool)
idx = np.argsort(D, axis=1)[:, :31]
for p in range(len(X)):
    L[p, idx[p]] = True

# 12x12 Jacobian of the point model, treated as one dense block (exact derivative, LU)
J = np.ones((12, 12), bool)

fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.45), gridspec_kw=dict(width_ratios=[1, 1, 0.62]))
ax[0].spy(K, markersize=0.15, color="#1d4e89")
ax[0].set_title(f"(a) global stiffness matrix\nonly $\\mathbf{{u}}$: {3*N}$\\times${3*N}, solved", fontsize=8, loc="left")
ax[1].spy(L, markersize=0.15, color="#b4530a")
ax[1].set_title(f"(b) field update at {len(X)} Gauss points\n$\\phi_{{3D}}, c, a$: 30 neighbours, explicit", fontsize=8, loc="left")
ax[2].spy(J, markersize=4.2, color="#3e4c59")
ax[2].set_title("(c) point model, one point\n12$\\times$12 Newton, local", fontsize=8, loc="left")
lab = [r"$\phi_{1..5}$", r"$\phi_0$", r"$\psi_{1..5}$", r"$\gamma$"]
ax[2].set_xticks([2, 5, 8, 11]); ax[2].set_xticklabels(lab, fontsize=6.5)
ax[2].set_yticks([2, 5, 8, 11]); ax[2].set_yticklabels(lab, fontsize=6.5)
for a in ax[:2]:
    a.set_xticks([]); a.set_yticks([])
fig.tight_layout()
A = HERE.parents[1] / "assets"
fig.savefig(A / "fig_matrix_structure.pdf"); fig.savefig(A / "fig_matrix_structure.png", dpi=400)
print("wrote", A / "fig_matrix_structure.png")
