"""anim_maze_rules.py -- video of the maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11) with the agreed coupling and the
antibiotic from T* = 20, for the two rules for newly reached points (A-73, A-74): left, the state of the neighbours
(the biofilm recovers); right, the start state (it does not). Top: living volume fraction sum_i phibar_i where the
biofilm phi_3D >= 0.05; bottom: antibiotic a. Python finite volumes, 4 um cells (maze_coupled_2d.py --snaps 120;
the Abaqus element gives the same, A-74); walls grey. A demonstration, not a thesis figure.

    python papers/coupling_paper/anim_maze_rules.py -> assets/anim_maze_rules.mp4, assets/anim_maze_rules_end.png
"""
import sys
from pathlib import Path

import imageio.v2 as imageio
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402

figstyle.apply()
Z = {"inherit": np.load(ROOT / "figdata" / "maze_h4_tend60_ab_k24fig_snaps.npz"),
     "fresh": np.load(ROOT / "figdata" / "maze_h4_tend60_ab_k24fig_fresh_snaps.npz")}
z0 = Z["inherit"]
X, Y, OPEN, NUT = z0["X"], z0["Y"], z0["open"].astype(bool), z0["nut"].astype(bool)   # um, indexing ij
ext = [X.min(), X.max(), Y.min(), Y.max()]
times = sorted({float(k.split("_", 1)[1]) for k in z0.files if k.startswith("phi_")})
TITLE = {"inherit": "new points take their neighbours' state", "fresh": "new points take the start state"}
CAPTION = (r"Maze of Klempt et al. 2024 (Sect. 4.4, Fig. 11) with the coupled model: Klempt et al. 2024 Eqs. 34$-$36 with the "
           r"point model of Klempt et al. 2026" "\n"
           r"(Fritsch et al. 2025, Case I) at every cell; nutrient $c = 1$ in the centre, colony at the entrance, antibiotic "
           r"$a = 1$ at the entrance from $T^* = 20$ ($\alpha^*_0 = 100$)." "\n"
           r"Table 6 of Klempt et al. 2024: $\beta = 10\,\mu$m$^2/T^*$, $r = 200\,\mu$m$/T^*$, $k_\alpha = 10^{-6}/T^*$; "
           r"4 $\mu$m cells. Walls grey, nutrient source outlined in blue. The two columns differ only in the rule for newly reached points.")

fig = plt.figure(figsize=(8.0, 8.4), dpi=160)  # 1280 x 1344 px
axs = [[fig.add_axes([0.03 + 0.48 * j, 0.53 - 0.405 * i, 0.40, 0.37]) for j in range(2)] for i in range(2)]
cax = [fig.add_axes([0.89, 0.53 - 0.405 * i, 0.018, 0.37]) for i in range(2)]
out = ROOT / "assets" / "anim_maze_rules.mp4"
movie = imageio.get_writer(str(out), fps=10, quality=9)
for T in times:
    for t in list(fig.texts):
        t.remove()
    for j, rule in enumerate(("inherit", "fresh")):
        z = Z[rule]
        phi, L, a = z[f"phi_{T}"], z[f"L_{T}"], z[f"a_{T}"]
        for i, (fld, cmap, mask) in enumerate(((L, "Greens", OPEN & (phi >= 0.05)), (a, "Reds", OPEN))):
            ax = axs[i][j]
            ax.clear()
            img = np.where(mask, fld, np.nan)
            m = ax.imshow(img.T, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
            ax.imshow(np.where(OPEN, np.nan, 1.0).T, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=4,
                      interpolation="nearest")
            ax.contour(X[:, 0] if X.ndim == 2 else X, Y[0, :] if Y.ndim == 2 else Y, NUT.T.astype(float), [0.5],
                       colors="tab:blue", linewidths=0.6)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(f"({'ab'[j]}) {TITLE[rule]}", fontsize=10, loc="left")
            if j == 1:
                c = fig.colorbar(m, cax=cax[i])
                c.set_label(r"living volume fraction $\Sigma_i\,\bar\phi_i$" if i == 0 else "antibiotic $a$", fontsize=9)
                c.ax.tick_params(labelsize=8)
    fig.text(0.97, 0.985, rf"$T^* = {T:.1f}$", ha="right", va="top", fontsize=12)
    fig.text(0.03, 0.105, CAPTION, ha="left", va="top", fontsize=7.6, linespacing=1.45)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    for _ in range(15 if T == times[-1] else 1):
        movie.append_data(img)
    if T == times[-1]:
        imageio.imwrite(out.with_name(out.stem + "_end.png"), img)
movie.close()
print("wrote", out)
