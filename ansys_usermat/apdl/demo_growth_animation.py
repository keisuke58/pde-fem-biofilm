"""demo_growth_animation.py -- a short animation for the colloquium: a biofilm seed on the bottom face of the 2 mm
cube spreads, grows and builds up stress. Illustration only, not a thesis result.

Model (simplified, parameters chosen for the picture, not calibrated):
  - phase phi on the 16^3 element grid: d phi/dt = beta lap(phi) + r phi (1 - phi), no-flux faces,
    beta = 0.02 mm^2/T* (the value of the ANSYS runs), r = 3 /T*;
  - growth alpha - 1 = k_alpha * int phi dt (Klempt 2024 Eq. 36 with c = 1), k_alpha = 0.02 /T*;
  - small-strain hex8 solve per frame with the eigenstrain (alpha - 1) I, stiffness E (phi^2 + f), f = 1e-3,
    nu = 0.49 (B-bar, element_matrices of mesh_study_seed.py), bottom face z = 0 clamped (the tooth surface).
Panels: the vertical mid-plane y = 1 mm of phi, alpha - 1 and von Mises / E; the outline shows the displacement,
magnified.

    python ansys_usermat/apdl/demo_growth_animation.py   -> assets/demo_growth.gif (and the last frame as .png)
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, str(HERE))
from mesh_study_seed import element_matrices  # noqa: E402

N, L = 16, 2.0
h = L / N
BETA, R, KA, FLOOR, NU = 0.02, 3.0, 0.02, 1e-3, 0.49
T_END, DT, NFRAMES = 3.0, 0.005, 60
MAG = 8.0                                               # displacement magnification in the outline

# ---- phase field on the element grid (finite volumes, explicit)
ax = h * (np.arange(N) + 0.5)
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
phi = np.where((np.abs(X - 1) < 0.25) & (np.abs(Y - 1) < 0.25) & (Z < 0.25), 1.0, 0.0)
alpha = np.zeros_like(phi)                              # alpha - 1


def lap(f):
    p = np.pad(f, 1, mode="edge")                       # no-flux faces
    return (p[2:, 1:-1, 1:-1] + p[:-2, 1:-1, 1:-1] + p[1:-1, 2:, 1:-1] + p[1:-1, :-2, 1:-1]
            + p[1:-1, 1:-1, 2:] + p[1:-1, 1:-1, :-2] - 6 * f) / h ** 2


# ---- hex8 mechanics
Ke0, fe0, B0, D, m = element_matrices(h, NU)
idx = np.arange((N + 1) ** 3).reshape(N + 1, N + 1, N + 1)
off = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
conn = np.stack([idx[i:i + N, j:j + N, k:k + N].ravel() for i, j, k in off], 1)
dofs = np.concatenate([3 * conn[:, [a]] + np.arange(3) for a in range(8)], 1)
rows = np.repeat(dofs, 24, 1).ravel(); cols = np.tile(dofs, (1, 24)).ravel()
ndof = 3 * (N + 1) ** 3
fix = (3 * idx[:, :, 0].ravel()[:, None] + np.arange(3)).ravel()
free = np.setdiff1d(np.arange(ndof), fix)


def mechanics(phi, alpha):
    Ee = (phi.ravel() ** 2 + FLOOR)                     # E = 1: stresses are von Mises / E
    eps0 = alpha.ravel()
    K = sp.csr_matrix(((Ke0[None] * Ee[:, None, None]).ravel(), (rows, cols)), shape=(ndof, ndof))
    F = np.zeros(ndof)
    np.add.at(F, dofs.ravel(), (fe0[None] * (Ee * eps0)[:, None]).ravel())
    u = np.zeros(ndof)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free])
    eps = (B0 @ u[dofs].T).T - eps0[:, None] * m
    sig = Ee[:, None] * (eps @ D.T)
    s = sig[:, :3] - sig[:, :3].mean(1, keepdims=True)
    vm = np.sqrt(1.5 * ((s ** 2).sum(1) + 2 * (sig[:, 3:] ** 2).sum(1)))
    return vm.reshape(N, N, N), u.reshape(N + 1, N + 1, N + 1, 3)


# ---- march and store frames
steps = int(round(T_END / DT))
every = steps // NFRAMES
frames = []
for n in range(steps + 1):
    if n % every == 0:
        vm, u = mechanics(phi, alpha)
        j = N // 2
        frames.append((n * DT, phi[:, j, :].copy(), alpha[:, j, :].copy(), vm[:, j, :].copy(),
                       u[:, j, :, :].copy()))
    phi = np.clip(phi + DT * (BETA * lap(phi) + R * phi * (1 - phi)), 0, 1)
    alpha = alpha + DT * KA * phi
vmax_vm = max(f[3].max() for f in frames)
amax = max(f[2].max() for f in frames)

# ---- draw
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9))
edges = np.linspace(0, L, N + 1)
xn = np.linspace(0, L, N + 1)
XN, ZN = np.meshgrid(xn, xn, indexing="ij")
specs = [(1, "phase $\\phi$", "Greens", 1.0), (2, "growth $\\alpha-1$", "Oranges", amax),
         (3, "von Mises / $E$", "viridis", vmax_vm)]
ims, lines = [], []
for a, (k, lab, cm, vmax) in zip(axs, specs):
    im = a.pcolormesh(edges, edges, frames[0][k].T, cmap=cm, vmin=0, vmax=vmax, shading="flat")
    cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7.5)
    if k == 2:
        cb.formatter.set_powerlimits((-2, 2)); cb.update_ticks()
    a.set_title(lab, fontsize=10, pad=3)
    a.set_aspect("equal"); a.set_xlim(-0.05, L + 0.25); a.set_ylim(-0.12, L + 0.25)
    a.set_xticks([0, 1, 2]); a.set_yticks([0, 1, 2])
    a.tick_params(labelsize=8); a.grid(False)
    a.axhline(0, color="0.3", lw=2.5, solid_capstyle="butt")
    ln, = a.plot([], [], color="k", lw=0.8)
    ims.append(im); lines.append(ln)
axs[0].set_ylabel("height $z$ (mm)")
for a in axs:
    a.set_xlabel("$x$ (mm)")
txt = fig.text(0.5, 0.015, "", ha="center", fontsize=9)
fig.text(0.995, 0.015, "illustration, parameters not calibrated", ha="right", fontsize=7, color="0.4")
fig.tight_layout(rect=(0, 0.05, 1, 1))


def outline(u):
    xs = XN + MAG * u[:, :, 0]; zs = ZN + MAG * u[:, :, 2]
    bx = np.concatenate([xs[:, 0], xs[-1, :], xs[::-1, -1], xs[0, ::-1]])
    bz = np.concatenate([zs[:, 0], zs[-1, :], zs[::-1, -1], zs[0, ::-1]])
    return bx, bz


def draw(i):
    t, *fields, u = frames[i]
    for im, f in zip(ims, fields):
        im.set_array(f.T.ravel())
    bx, bz = outline(u)
    for ln in lines:
        ln.set_data(bx, bz)
    txt.set_text(f"$t$ = {t:.2f} $T^*$   (bottom face fixed, outline displacement x{MAG:g})")
    return ims + lines + [txt]


anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
out = ROOT / "assets" / "demo_growth"
anim.save(out.with_suffix(".gif"), writer=PillowWriter(fps=12), dpi=110)
draw(len(frames) - 1)
fig.savefig(out.with_suffix(".png"), dpi=170)
print("wrote", out.with_suffix(".gif"), len(frames), "frames; max vM/E", round(vmax_vm, 4), "max alpha-1",
      round(amax, 4))
