"""validation_volume.py -- the model's biofilm amount V(t) = int phi_3D dV, divided by its value at day 1, at the data
days of Heine et al. 2025 (Front. Oral Health 6:1649419, Fig. 2A; days 1, 3, 6, 10, 15, 21), from the layer runs of
coupled_living_consumption.py --example layer --eco-file (papers/coupling_paper/VALIDATION_PLAN.ja.md). T* = 1 spans
day 1 to day 21, so day d is T* = (d - 1)/20.

    python ansys_usermat/apdl/validation_volume.py
        -> ansys_usermat/apdl/results/2026-10-validation/volume_rel_model.json (model values only)
    python ansys_usermat/apdl/validation_volume.py --measured heine_volume_rel.json --out fig.pdf
        also draws the comparison (model curves, measured medians with their 90 % bootstrap band). The measured
        values are not part of this repository (the plan: no raw data here); the figure goes where --out points.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "ansys_usermat"))
RES = HERE / "results" / "2026-10-validation"
DAYS = [1, 3, 6, 10, 15, 21]
COND = ["CH", "CS", "DH", "DS"]


def opt(name, default=None):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def load_runs():
    out = {}
    for f in sorted(RES.glob("coupled_living_consumption_*_layer.json")):
        d = json.loads(f.read_text())
        h = d["runs"]["layer"]
        m = re.search(r"_ecocal5_(\w\w)", f.name)
        if not m:
            continue
        law = "b" if "_growlivE" in f.name else "a"
        nf = int(re.search(r"_nf(\d+)", f.name).group(1)) if "_nf" in f.name else 32
        h0 = float(re.search(r"_h0([\d.]+)", f.name).group(1)) if "_h0" in f.name else 0.25
        t = np.asarray(h["t"]); V = np.asarray(h["phi_total_mm3"]); L = np.asarray(h["living_mean"])
        Td = [(dd - 1) / 20.0 for dd in DAYS]
        out[(m.group(1), law, nf, h0)] = dict(
            file=f.name, condition=m.group(1), law=law, nf=nf, layer_h0_frac=h0, setting=d["setting"],
            days=DAYS, T_star=Td, V_rel=np.interp(Td, t, V / V[0]).tolist(),
            living_mean=np.interp(Td, t, L).tolist(), V_day1_mm3=float(V[0]),
            V_rel_end=float(V[-1] / V[0]), t=t.tolist(), V_rel_t=(V / V[0]).tolist(), living_mean_t=L.tolist())
    return out


def main():
    runs = load_runs()
    LAW = {"a": "published law (Klempt 2024 Eq. 34-36)",
           "b": "variant A, --growth-living-energy (-k_alpha alpha phi sum phibar_i)"}
    res = {"description": "V(t)/V(day 1), V = int phi_3D dV, and the phi-weighted mean living fraction sum_i phibar_i "
                          "of the layer runs at the data days of Heine et al. 2025 (1, 3, 6, 10, 15, 21); "
                          "T* = (day - 1)/20. Model values only; the measured values stay outside the repository.",
           "laws": LAW, "runs": []}
    for k in sorted(runs):
        r = runs[k]
        res["runs"].append({kk: r[kk] for kk in ("condition", "law", "nf", "layer_h0_frac", "file", "days", "T_star",
                                                  "V_rel", "living_mean", "V_day1_mm3", "V_rel_end")})
        print(f"{r['condition']} law {r['law']} {r['nf']}^3 h0 {r['layer_h0_frac']:g}: V/V1 at days "
              + ", ".join(f"{dd}: {v:.3f}" for dd, v in zip(DAYS, r["V_rel"]))
              + " | living " + ", ".join(f"{v:.2f}" for v in r["living_mean"]))
    (RES / "volume_rel_model.json").write_text(json.dumps(res, indent=1) + "\n")
    print("wrote", RES / "volume_rel_model.json")
    meas = opt("--measured")
    if not meas:
        return
    M = json.loads(Path(meas).read_text())
    import figstyle
    figstyle.apply(size=8)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nf_main = int(opt("--nf", 32))
    fig, axs = plt.subplots(2, 2, figsize=(6.3, 4.6), sharex=True)
    for ax, C in zip(axs.ravel(), COND):
        days = [int(d) for d in M[C]]
        rel = [M[C][str(d)]["rel"] for d in days]
        lo = [M[C][str(d)]["lo"] for d in days]; hi = [M[C][str(d)]["hi"] for d in days]
        ax.errorbar(days, rel, yerr=[np.subtract(rel, lo), np.subtract(hi, rel)], fmt="o", ms=3.5, color="k",
                    capsize=2, lw=0.8, label="measured (median, 90 % band)", zorder=5)
        for law, col, ls in (("a", "C0", "-"), ("b", "C3", "--")):
            k = (C, law, nf_main, 0.25)
            if k in runs:
                r = runs[k]
                ax.plot(1 + 20 * np.asarray(r["t"]), r["V_rel_t"], color=col, ls=ls, lw=1.2,
                        label={"a": "model, published law", "b": "model, variant A"}[law])
            k2 = (C, law, nf_main, 0.125)
            if k2 in runs:
                r = runs[k2]
                ax.plot(1 + 20 * np.asarray(r["t"]), r["V_rel_t"], color=col, ls=":", lw=0.9,
                        label={"a": "published law, $h_0$ = 0.125 box", "b": "variant A, $h_0$ = 0.125 box"}[law])
        ax.set_title(C, fontsize=8, loc="left")
        ax.set_xticks(DAYS); ax.grid(False)
        ax.set_yscale("log")
    for ax in axs[1]:
        ax.set_xlabel("day")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"$V/V_{\mathrm{day\,1}}$")
    axs[0, 0].legend(fontsize=6, frameon=False, loc="upper left")
    fig.tight_layout()
    out = Path(opt("--out", "validation_volume.pdf"))
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=200, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
