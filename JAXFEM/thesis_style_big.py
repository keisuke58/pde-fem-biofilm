"""thesis_style_big.py — slide / poster variant of thesis_style.

Same usetex preamble and palette as thesis_style.py, but font sizes scaled
for A0 posters and Beamer presentations where text must be legible at 3+ m.

Usage in a figure/slide script:
    from thesis_style_big import use as thesis_use, clean_ax, PALETTE
    figsize = thesis_use(width_frac=1.0, aspect=0.56)
    fig, ax = plt.subplots(figsize=figsize)
    ...

Differences from thesis_style.py:
  * font.size   9  → 18  (labels, legends, text)
  * axes title  9  → 20
  * tick labels 8  → 16
  * lines       1.2 → 2.0
  * axes lines  0.6 → 1.2
  * savefig dpi 300 → 150  (slides don't need print resolution)
  * TEXTWIDTH_IN based on Beamer 16:9 frame (\textwidth ≈ 12.8 cm = 5.04 in)

Width fractions (Beamer metropolis 16:9, inner margins ~12.8 cm usable):
  width_frac=1.0  → full-width column figure
  width_frac=0.48 → two-column layout (two figures side by side)
"""
import matplotlib as mpl
from thesis_style import _PREAMBLE, PALETTE, SP_COLORS, SP_SHORT, clean_ax, tight_ylim

BEAMER_TEXTWIDTH_IN = 5.04   # 12.8 cm usable width in metropolis 16:9

def use(width_frac=1.0, aspect=0.62):
    """Apply big-font rcParams and return figsize (w, h) in inches."""
    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": [],
        "text.latex.preamble": _PREAMBLE,
        "font.size": 18,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "figure.titlesize": 20,
        "lines.linewidth": 2.0,
        "axes.linewidth": 1.2,
        "grid.linewidth": 0.8,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    w = BEAMER_TEXTWIDTH_IN * width_frac
    return (w, w * aspect)
