"""klempt_table3_fig.py -- Klempt et al. 2024 Table 3 (test case 4.1) from the
Abaqus run of make_klempt_inp.py --case fig4_corner --blend 0.5 --first-order
--E 10 --s-every 25, next to the paper's panels: phi, c and the hydrostatic
stress p on the diagonal cut x = z (nutrient corner at the top right), rows
5 / 25 / 45 / 100 %.

    python abaqus_composition/klempt_table3_fig.py JOB.dat [OUT.png]  (default assets/fig_abaqus_klempt2024_table3.png)

Element centroid values on the 20^3 mesh. phi and c are printed every 10
increments, so at 5 / 25 / 45 % they are the mean of the two neighbouring
prints (+-0.01); p is printed at these times. p in the unit of E, read as MPa
(E = 10, KLEMPT2024_REPRODUCTION.md sec. 13). Colours: the paper's banded
legends (figstyle.klempt_cmap); p mapped on +3e-4 ... -6.5e-4 MPa.
The paper's panels are cut from the PDF (CC BY 4.0, see THIRD_PARTY.md).
"""
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
sys.path[:0] = [str(HERE), str(ROOT / "ansys_usermat")]
import figstyle  # noqa: E402
from table3_pressure import tables  # noqa: E402

OUT = ROOT / "assets" / "fig_abaqus_klempt2024_table3.png"
PDF = ROOT / "references" / "Klempt2024_Hamilton_biofilm_growth_BMMB.pdf"
TIMES = [0.05, 0.25, 0.45, 1.00]
PAPER_COLS = {"phi": (957, 1287), "c": (1402, 1733), "p": (1847, 2178)}
PAPER_ROWS = [(291, 524), (536, 770), (782, 1015), (1028, 1262)]
P_HI, P_LO = 3e-4, -6.5e-4
N = 20


def grid(tab, name):
    u = np.zeros((N, N, N))
    i = tab["head"].index(name)
    for e, row in tab["rows"].items():
        q = e - 1
        u[q % N, (q // N) % N, q // N // N] = row[i]
    return u


def fields(dat):
    T = tables(dat)
    get = lambda t, name: next(x for x in T[round(t, 6)] if name in x["head"])  # noqa: E731
    out = {}
    for t in TIMES:
        ts = [t] if round(t, 6) in T and any("TEMP" in x["head"] for x in T[round(t, 6)]) else [t - 0.01, t + 0.01]
        ts = [x for x in ts if round(x, 6) in T] or [t - 0.01]
        phi = np.mean([grid(get(x, "TEMP"), "TEMP") for x in ts], axis=0)
        c = np.mean([grid(get(x, "TEMP"), "SDV51") for x in ts], axis=0)
        S = get(t, "S11")
        p = (grid(S, "S11") + grid(S, "S22") + grid(S, "S33")) / 3
        out[t] = (phi, c, p)
    return out


def diagonal_cut(u):
    return np.array([[u[i, j, i] for i in range(N)] for j in range(N)])


def paper_panels():
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "11", "-l", "11", "-r", "300", "-png", str(PDF),
                        str(Path(tmp) / "p")], check=True)
        im = Image.open(next(Path(tmp).glob("p*.png"))).convert("RGB")
        return {(k, r): np.asarray(im.crop((x0, y0, x1, y1)))
                for k, (x0, x1) in PAPER_COLS.items() for r, (y0, y1) in enumerate(PAPER_ROWS)}


def main(dat):
    figstyle.apply(size=11)
    cmap, norm = figstyle.klempt_cmap()
    F = fields(dat)
    P = paper_panels()
    L = 20.0
    fig, ax = plt.subplots(len(TIMES), 6, figsize=(17, 9.2))
    heads = [r"$\phi$, Abaqus", r"$\phi$, paper", r"$c$, Abaqus", r"$c$, paper",
             r"$p$, Abaqus", r"$p$, paper"]
    for r, t in enumerate(TIMES):
        phi, c, p = F[t]
        pn = (p - P_LO) / (P_HI - P_LO)
        cells = (("sim", phi), ("pap", "phi"), ("sim", c), ("pap", "c"), ("sim", pn), ("pap", "p"))
        for col, (kind, data) in enumerate(cells):
            a = ax[r, col]
            if kind == "sim":
                a.imshow(diagonal_cut(data), origin="lower", extent=(0, L * np.sqrt(2), 0, L),
                         cmap=cmap, norm=norm, interpolation="nearest", aspect="auto")
            else:
                a.imshow(P[(data, r)], extent=(0, L * np.sqrt(2), 0, L), aspect="auto")
            a.set_xticks([]); a.set_yticks([]); a.grid(False)
            if r == 0:
                a.set_title(heads[col])
        ax[r, 0].set_ylabel(f"{round(t * 100)} %")
        ax[r, 4].text(0.02, 0.02, f"{p.min():.1e} ... {p.max():.1e}", transform=ax[r, 4].transAxes,
                      fontsize=8, color="k", bbox=dict(fc="w", ec="none", alpha=0.7))
    fig.suptitle("Klempt et al. 2024, test case 4.1 (Table 3): diagonal cut, nutrient corner at the top right. "
                 "Abaqus 20$^3$: growth on every face (w = 0.5), consumption $g\\phi c$,\nnutrient in a 2 um corner "
                 "block, E = 10 read as MPa (assumptions of this work); p colours on +3e-4 ... -6.5e-4 MPa as the "
                 "paper's legend. Paper: Table 3 (CC BY 4.0)", fontsize=11)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        OUT = Path(sys.argv[2])
    main(sys.argv[1])
