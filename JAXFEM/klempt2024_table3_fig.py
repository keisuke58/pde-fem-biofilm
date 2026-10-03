#!/usr/bin/env python3
"""Klempt et al. 2024 Table 3 (test case 4.1, growth from the centre towards
the nutrient corner), redrawn from this repository's reproduction next to the
paper's own panels.

Field: klempt2024_timescale.py's fit (Table 2 unchanged, nutrient on the edge
line, time scale 1.5, growth on both faces, consumption g phi c). Cut: the
diagonal plane through the middle that the paper uses, with the nutrient
corner at the top right (aspect sqrt(2) : 1). Colours: the paper's banded
legend (figstyle.klempt_cmap). The paper's panels are cut from the PDF
(CC BY 4.0, see THIRD_PARTY.md).

    python JAXFEM/klempt2024_table3_fig.py   -> assets/fig_klempt2024_table3.png
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "ansys_usermat"))
import figstyle  # noqa: E402
import klempt2024_case1_bc as B  # noqa: E402
import klempt2024_quantitative as K  # noqa: E402

OUT = ROOT / "assets" / "fig_klempt2024_table3.png"
PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
TIMES = [0.05, 0.25, 0.45, 1.00]
S = 1.5
# Table 3 on page 11 at 300 dpi: phi and c columns, rows 5 / 25 / 45 / 100 %
PAPER_COLS = {"phi": (957, 1287), "c": (1402, 1733)}
PAPER_ROWS = [(291, 524), (536, 770), (782, 1015), (1028, 1262)]


def fields():
    K.R, K.BETA, K.K_A = 100 * S, 2 * S, 1e-3 * S
    phi, _, _ = K.setup("fig4_edge")
    mask = B.nutrient_mask("edge")
    out, t = {}, 0.0
    for t_next in TIMES:
        _, phi, c = B.run_setup(phi, mask, 1e8, "ic", "first_order", "abs",
                                t_end=round(t_next - t, 6))
        out[t_next] = (phi.copy(), c.copy())
        t = t_next
    return out


def diagonal_cut(u):
    """The plane x = z (contains y): horizontal = the x = z diagonal, vertical
    = y. The nutrient edge x = y = L meets it at (L, L, L): the top right."""
    n = u.shape[0]
    return np.array([[u[i, j, i] for i in range(n)] for j in range(n)])


def paper_panels():
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "11", "-l", "11", "-r", "300", "-png", str(PDF),
                        str(Path(tmp) / "p")], check=True)
        im = Image.open(next(Path(tmp).glob("p*.png"))).convert("RGB")
        return {(k, r): np.asarray(im.crop((x0, y0, x1, y1)))
                for k, (x0, x1) in PAPER_COLS.items() for r, (y0, y1) in enumerate(PAPER_ROWS)}


def main():
    figstyle.apply(size=11)
    cmap, norm = figstyle.klempt_cmap()
    F = fields()
    P = paper_panels()
    L = K.L
    fig, ax = plt.subplots(len(TIMES), 4, figsize=(12.5, 9.2))
    heads = [r"$\phi$, this work", r"$\phi$, paper", r"$c$, this work", r"$c$, paper"]
    for r, t in enumerate(TIMES):
        phi, c = F[t]
        for col, (kind, data) in enumerate((("sim", phi), ("pap", "phi"), ("sim", c), ("pap", "c"))):
            a = ax[r, col]
            if kind == "sim":
                a.imshow(diagonal_cut(data), origin="lower", extent=(0, L * np.sqrt(2), 0, L),
                         cmap=cmap, norm=norm, interpolation="bilinear", aspect="auto")
                figstyle.element_grid(a, np.arange(0, 20.5) * np.sqrt(2), np.arange(0, 20.5))
            else:
                a.imshow(P[(data, r)], extent=(0, L * np.sqrt(2), 0, L), aspect="auto")
            a.set_xticks([]); a.set_yticks([]); a.grid(False)
            for sp in a.spines.values():
                sp.set_visible(True)
            if r == 0:
                a.set_title(heads[col])
        ax[r, 0].set_ylabel(f"{round(t * 100)} %")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, label=r"$\phi$, $c$  [-]")
    fig.suptitle("Klempt et al. 2024, test case 4.1 (Table 3): diagonal cut, nutrient corner at the "
                 "top right.\nThis work: Python reproduction (Table 2, time scale 1.5, both-face growth, "
                 r"consumption $g\phi c$); paper: Table 3 (CC BY 4.0)", fontsize=11.5)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
