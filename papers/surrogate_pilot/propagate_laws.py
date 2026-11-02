#!/usr/bin/env python3
"""The same uncertainty as propagate.py carried through the two growth laws, both on 64^3 at the same 16 points:
the published law (runs f) and the growth of the living biomass (runs e, --growth-living-energy, paper Sect. 2.5).
Keio continuation, not the thesis. -> assets/fig_surrogate_laws.png/.pdf, propagation_laws.json

    python papers/surrogate_pilot/propagate_laws.py
"""
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
sys.path.insert(0, str(HERE))
import figstyle  # noqa: E402
from fit_pilot import LO, HI, gp, load  # noqa: E402

warnings.filterwarnings("ignore")
N = 10_000


def main():
    figstyle.apply()
    rng = np.random.default_rng(0)
    a0 = np.clip(np.exp(np.log(100.0) + rng.normal(0.0, np.log(2.0) / 2, N)), 10.0, 300.0)
    gc = rng.uniform(0.7, 1.3, N)
    Xs = (np.c_[np.log10(a0), gc] - LO) / (HI - LO)
    out, samples = {}, {}
    for key, tag, lab in (("published", "f", "published law (Eq. 36)"), ("living", "e", "growth of the living biomass")):
        X, _, V, _ = load(tag)
        g = gp().fit(X, V); mu, sd = g.predict(Xs, return_std=True)
        Vs = mu + sd * rng.normal(size=N)
        samples[key] = (Vs, lab)
        p = [float(np.percentile(Vs, q)) for q in (5, 50, 95)]
        out[key] = {"vm_1e-4Pa_5_50_95": p, "spread_90_percent_of_median": (p[2] - p[0]) / p[1],
                    "runs_min_max": [float(V.min()), float(V.max())]}
    (HERE / "propagation_laws.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    fig, ax = plt.subplots(figsize=(3.8, 2.4), constrained_layout=True)
    for (Vs, lab), col in zip(samples.values(), ("#1f4e8c", "#c92a2a")):
        ax.hist(Vs, bins=60, color=col, alpha=0.55, label=lab)
    ax.set_xlabel(r"seed $\sigma_\mathrm{vM}$ [$10^{-4}$ Pa], $64^3$"); ax.set_ylabel("samples")
    ax.legend(frameon=False, fontsize=7.5)
    o = ROOT / "assets" / "fig_surrogate_laws.png"
    fig.savefig(o, dpi=600); fig.savefig(o.with_suffix(".pdf"))
    print("wrote", o)


if __name__ == "__main__":
    main()
