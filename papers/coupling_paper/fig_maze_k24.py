"""fig_maze_k24.py -- the maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11, read off the figure: abaqus_living/maze_k24.txt)
with the agreed coupling, ansys_usermat/apdl/maze_coupled_2d.py --maze-file ... --tend 60 [--ab]: phi_3D at four times,
without (top) and with the antibiotic (bottom, a = 1 at the entrance from T* = 20). Walls white, nutrient source
outlined in red.

    python papers/coupling_paper/fig_maze_k24.py [--h 2] -> papers/coupling_paper/fig_maze_k24_h<h>.pdf
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply()
H = float(sys.argv[sys.argv.index("--h") + 1]) if "--h" in sys.argv else 2.0
TIMES = (12.0, 24.0, 36.0, 60.0)
fig, axs = plt.subplots(2, 4, figsize=(6.3, 3.45), constrained_layout=True)
for row, ab in enumerate(("", "_ab")):
    z = np.load(ROOT / "figdata" / f"maze_h{H:g}_tend60{ab}_k24fig_transc.npz")
    op, nut = z["open"], z["nut"]
    for col, t in enumerate(TIMES):
        ax = axs[row, col]
        ph = np.ma.masked_where(~op, z[f"phi_{t:.1f}"])
        # rows of the mask file = rows of the paper's figure (entrance at the top)
        im = ax.imshow(ph, origin="upper", cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
        ax.contour(nut.astype(float), levels=[0.5], colors="r", linewidths=0.6)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$T^*={t:g}$", fontsize=8)
        if col == 0:
            ax.set_ylabel("no antibiotic" if not ab else "antibiotic", fontsize=8)
    axs[row, 0].text(-0.08, 1.05, "(a)" if row == 0 else "(b)", transform=axs[row, 0].transAxes, fontsize=9)
cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01)
cb.set_label(r"$\phi_{3D}$", fontsize=8)
out = HERE / f"fig_maze_k24_h{H:g}.pdf"
fig.savefig(out)
print("wrote", out)
