"""maze_from_figure.py -- the maze of Klempt et al. 2024, Sect. 4.4, Fig. 11 (CC BY 4.0), read off the figure by colour
into a cell mask (A-49): grey = the path (discretised body), white = walls (not discretised, 10 um) and outside, red = the
initial biofilm (phi = 1) at the entrance, green = the nutrient source (c = 1) in the centre. Diameter 560 um; one cell
layer of 4 um (the paper: semi-2D, one element thick, 4 um).

    pdftoppm -f 18 -l 18 -r 600 -png references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf <dir>/p
    python abaqus_living/maze_from_figure.py <dir>/p-18.png [--cell 4] [--out abaqus_living/maze_k24]

Writes <out>.npz (mask: 0 wall/outside, 1 path, 2 initial biofilm, 3 nutrient source; cell size in um) and <out>.png.
"""
import argparse
import sys

import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("png")
ap.add_argument("--cell", type=float, default=4.0, help="cell size in um")
ap.add_argument("--out", default="abaqus_living/maze_k24")
a = ap.parse_args()

im = np.asarray(Image.open(a.png).convert("RGB")).astype(int)
R, G, B = im[..., 0], im[..., 1], im[..., 2]
grey = (abs(R - 120) < 30) & (abs(G - 135) < 30) & (abs(B - 155) < 30) & (B > R)
red = (R > 200) & (G < 80) & (B < 80)
green = (G > 200) & (R < 80) & (B < 80)
body = grey | red | green
# the maze: the largest blob of body pixels; its bounding box gives the 560 um diameter
ys, xs = np.nonzero(body)
# keep the figure area only (the maze is in the upper half of the page)
H, Wd = body.shape
sel = (ys < 0.45 * H) & (xs > 0.35 * Wd)
ys, xs = ys[sel], xs[sel]
y0, y1, x0, x1 = ys.min(), ys.max(), xs.min(), xs.max()
d_px = max(x1 - x0, y1 - y0) + 1
um_per_px = 560.0 / d_px
n = int(round(560.0 / a.cell))
mask = np.zeros((n, n), np.int8)
for i in range(n):          # i: row (y, from the top), j: column (x)
    for j in range(n):
        py0 = int(y0 + i * a.cell / um_per_px); py1 = max(py0 + 1, int(y0 + (i + 1) * a.cell / um_per_px))
        px0 = int(x0 + j * a.cell / um_per_px); px1 = max(px0 + 1, int(x0 + (j + 1) * a.cell / um_per_px))
        blk = (slice(py0, py1), slice(px0, px1))
        fr, fg, fb = red[blk].mean(), green[blk].mean(), grey[blk].mean()
        if fr > 0.5:
            mask[i, j] = 2
        elif fg > 0.5:
            mask[i, j] = 3
        elif fb + fr + fg > 0.5:
            mask[i, j] = 1
# the arrow of the figure crosses the green disk: every cell inside the disk is source
gi, gj = np.nonzero(mask == 3)
ci, cj = gi.mean(), gj.mean()
rad = np.sqrt(((gi - ci) ** 2 + (gj - cj) ** 2).max())
I, J = np.indices(mask.shape)
mask[((I - ci) ** 2 + (J - cj) ** 2 <= rad ** 2)] = 3
# no corner-only contact across a wall: two body cells that touch only at a corner (both other cells of the 2x2 block
# wall) would share a node in a finite-element mesh; the one with fewer body neighbours becomes wall
body = mask > 0
nfix = 0
for _ in range(10):
    changed = False
    for i in range(n - 1):
        for j in range(n - 1):
            b = body[i:i + 2, j:j + 2]
            for (p, q), (r, s) in ((((0, 0), (1, 1))), (((0, 1), (1, 0)))):
                if b[p, q] and b[r, s] and not b[p, s] and not b[r, q]:
                    cand = [(i + p, j + q), (i + r, j + s)]
                    nb = [sum(body[x + dx, y + dy] for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
                              if 0 <= x + dx < n and 0 <= y + dy < n) for x, y in cand]
                    x, y = cand[int(np.argmin(nb))]
                    if mask[x, y] in (2, 3):
                        x, y = cand[1 - int(np.argmin(nb))]
                    mask[x, y] = 0; body[x, y] = False; nfix += 1; changed = True
                    b = body[i:i + 2, j:j + 2]
    if not changed:
        break
print(f"corner-only contacts removed: {nfix}")
# connectivity through cell faces: the entrance and the source must lie in one component (the maze is solvable)
from collections import deque  # noqa: E402
comp = -np.ones(mask.shape, int); sizes = []
for s in zip(*np.nonzero(mask > 0)):
    if comp[s] >= 0:
        continue
    k = len(sizes); comp[s] = k; dq = deque([s]); cnt = 0
    while dq:
        x, y = dq.popleft(); cnt += 1
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            u, v = x + dx, y + dy
            if 0 <= u < n and 0 <= v < n and mask[u, v] > 0 and comp[u, v] < 0:
                comp[u, v] = k; dq.append((u, v))
    sizes.append(cnt)
print(f"components {len(sizes)} (largest {sorted(sizes)[-3:]}); entrance in {sorted(set(comp[mask == 2]))}, "
      f"source in {sorted(set(comp[mask == 3]))}")
np.savez_compressed(a.out + ".npz", mask=mask, cell_um=a.cell, diameter_um=560.0)
with open(a.out + ".txt", "w") as f:     # the committed form (npz and png are gitignored): one row per line, top first
    f.write(f"# maze of Klempt et al. 2024 Fig. 11, {n} x {n} cells of {a.cell:g} um; 0 wall/outside, 1 path, "
            f"2 initial biofilm, 3 nutrient source\n")
    f.writelines("".join(str(v) for v in row) + "\n" for row in mask)
col = np.array([[255, 255, 255], [122, 136, 156], [230, 0, 0], [0, 220, 0]], np.uint8)
Image.fromarray(col[mask]).resize((4 * n, 4 * n), Image.NEAREST).save(a.out + ".png")
print(f"{a.png}: maze {d_px} px = 560 um ({um_per_px:.3f} um/px); {n} x {n} cells of {a.cell} um: "
      f"path {int((mask == 1).sum())}, initial biofilm {int((mask == 2).sum())}, source {int((mask == 3).sum())}, "
      f"wall/outside {int((mask == 0).sum())}")
