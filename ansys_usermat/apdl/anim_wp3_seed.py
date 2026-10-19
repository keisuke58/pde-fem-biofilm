"""anim_wp3_seed.py -- animation of the verification problem of Chapter 4 (seed growth in the 2 mm cube, one species),
from the reference solution the element was checked against.

Field: wp3_reference.solve (finite volumes, RK4, beta = 0.02 mm^2/T*, k_alpha = 1e-3 /T*, penalty 5, the 32 seed
elements of the partner's 8^3 mesh, zero-flux walls) on 48^3, stored at 45 times up to T* = 1.1.
Stress: wp3_stress_reference.solve_fields (hex8, B-bar, E = (phi^2 + 1e-3) 10 Pa, nu = 0.49, eigenstrain (alpha - 1) I,
rigid-body constraints at the three deck corners) on 16^3 with the field averaged onto the elements; at T* = 1.1 the
seed von Mises on 16^3 is within 1 % of the extrapolated value (wp3_stress_reference.json).
Panels: the mid-plane z = 0 through the seed (mean of the two layers next to it): phi, alpha - 1, von Mises.

    python ansys_usermat/apdl/anim_wp3_seed.py   -> assets/anim_wp3_seed.gif, .png (last frame), .json (seed means)
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
import wp3_reference as W  # noqa: E402
import wp3_stress_reference as S  # noqa: E402

figstyle.apply(10)
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

NF, NM, NFRAMES = 48, 16, 88
cells8 = W.seed_cells_8(HERE / "results/2026-10-wp3_fix/wp3c_n8_dt00125.json")


def march(n, nframes):
    """wp3_reference.solve with the fields kept at nframes + 1 equally spaced times."""
    h = W.L / n
    dt = 0.15 * h * h / W.BETA
    per = int(np.ceil(W.T_END / nframes / dt))
    dt = W.T_END / (nframes * per)
    phi = W.seed_mask(n, cells8).astype(float)
    alpha = np.ones_like(phi)
    lap = np.empty_like(phi)
    out = [(0.0, phi.copy(), alpha - 1)]
    for f in range(nframes):
        for _ in range(per):
            k1p, k1a = W.rhs(phi, alpha, h, lap)
            k2p, k2a = W.rhs(phi + 0.5 * dt * k1p, alpha + 0.5 * dt * k1a, h, lap)
            k3p, k3a = W.rhs(phi + 0.5 * dt * k2p, alpha + 0.5 * dt * k2a, h, lap)
            k4p, k4a = W.rhs(phi + dt * k3p, alpha + dt * k3a, h, lap)
            phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
            alpha += dt / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        out.append(((f + 1) * per * dt, phi.copy(), alpha - 1))
    return out


def mid(field):
    k = field.shape[2] // 2
    return 0.5 * (field[:, :, k - 1] + field[:, :, k])


frames, seeds = [], []
seed16 = S.seed_mask(NM, cells8)
for t, phi, am1 in march(NF, NFRAMES):
    pm, am = W.block_mean(phi, NM), W.block_mean(am1, NM)
    vm, p = S.solve_fields(NM, pm, am)
    frames.append((t, mid(phi), mid(am1), mid(vm)))
    seeds.append({"t": round(t, 4), "seed alpha-1": float(am[seed16].mean()), "seed vM Pa": float(vm[seed16].mean()),
                  "seed p Pa": float(p[seed16].mean())})
    print(f"t={t:.3f}  seed vM {vm[seed16].mean():.3e} Pa", flush=True)
json.dump({"source": "anim_wp3_seed.py", "field grid": NF, "mech mesh": NM, "seed means": seeds},
          open(ROOT / "assets/anim_wp3_seed.json", "w"), indent=1)

amax = max(f[2].max() for f in frames)
vmax = max(f[3].max() for f in frames)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.75))
specs = [(1, "$\\phi$", "Greens", 1.0, NF, 1.0), (2, "$\\alpha-1$  ($10^{-3}$)", "Oranges", amax * 1e3, NF, 1e3),
         (3, "von Mises (mPa)", "viridis", vmax * 1e3, NM, 1e3)]
ims = []
for a, (k, lab, cm, vm_, n, sc) in zip(axs, specs):
    e = np.linspace(-1, 1, n + 1)
    im = a.pcolormesh(e, e, frames[0][k].T * sc, cmap=cm, vmin=0, vmax=vm_, shading="flat", rasterized=True)
    cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7.5)
    a.set_title(lab, fontsize=10, pad=3)
    a.set_aspect("equal"); a.set_xticks([-1, 0, 1]); a.set_yticks([-1, 0, 1])
    a.tick_params(labelsize=8); a.grid(False)
    a.set_xlabel("$x$ (mm)")
    ims.append((im, k, sc))
axs[0].set_ylabel("$y$ (mm)")
txt = fig.text(0.5, 0.02, "", ha="center", fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 1))


def draw(i):
    t = frames[i][0]
    for im, k, sc in ims:
        im.set_array((frames[i][k].T * sc).ravel())
    txt.set_text(f"$T^*$ = {t:.2f}   (mid-plane $z$ = 0 of the 2 mm cube; seed in the centre)")
    return [im for im, _, _ in ims] + [txt]


anim = FuncAnimation(fig, draw, frames=len(frames))
out = ROOT / "assets" / "anim_wp3_seed"
anim.save(out.with_suffix(".gif"), writer=PillowWriter(fps=16), dpi=220)
draw(len(frames) - 1)
fig.savefig(out.with_suffix(".png"), dpi=400); fig.savefig(out.with_suffix(".pdf"))
print("wrote", out.with_suffix(".gif"))
