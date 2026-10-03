"""Shared matplotlib style for thesis and slide figures: Times New Roman.

    import figstyle; figstyle.apply()

Times New Roman where it is installed (IKMHIWI03, Windows); elsewhere the
metric-compatible Liberation Serif / Nimbus Roman. Mathematics in STIX, which
matches Times.
"""
import matplotlib

SERIF = ["Times New Roman", "Liberation Serif", "Nimbus Roman", "TeX Gyre Termes",
         "STIXGeneral", "DejaVu Serif"]


def apply(size=11):
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": SERIF,
        "mathtext.fontset": "stix",
        "font.size": size,
        "axes.titlesize": size + 1,
        "axes.labelsize": size,
        "legend.fontsize": size - 1.5,
        "xtick.labelsize": size - 1,
        "ytick.labelsize": size - 1,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })


# Contour colours of Klempt et al. 2024 (Tables 3 and 4: ANSYS's default
# legend, 7 bands from blue = 0 to red = 1), read off the PDF. The legend is
# bright; the contour panels show the same bands shaded to about 2/3
# brightness, which is what the eye remembers of the paper.
KLEMPT_COLORS = ["#0000ff", "#00cbff", "#00ffcc", "#00ff00", "#cbff00", "#ffcb00", "#ff0000"]
KLEMPT_PANEL_COLORS = ["#0000a8", "#0177aa", "#00aa78", "#05a607", "#7bac05", "#ae7708", "#a80001"]


def klempt_cmap(vmin=0.0, vmax=1.0, shaded=True):
    """(cmap, norm) in the paper's banded colours; shaded=True gives the
    colours of its contour panels, False those of its legend."""
    import numpy as np
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cmap = ListedColormap(KLEMPT_PANEL_COLORS if shaded else KLEMPT_COLORS)
    return cmap, BoundaryNorm(np.linspace(vmin, vmax, cmap.N + 1), cmap.N)


def element_grid(ax, x_edges, y_edges, color="black", alpha=0.35, lw=0.3):
    """Thin element edges over a contour panel, as in the paper's figures."""
    for x in x_edges:
        ax.axvline(x, color=color, alpha=alpha, lw=lw)
    for y in y_edges:
        ax.axhline(y, color=color, alpha=alpha, lw=lw)
