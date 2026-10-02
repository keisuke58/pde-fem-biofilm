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
