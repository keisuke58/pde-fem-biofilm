#!/usr/bin/env python3
"""compare_cal5.py -- compare the calibrated five-species runs across the four
conditions (CH, CS, DS, DH) at the end of the run (T* = 1, Day 21).

    python ansys_usermat/apdl/compare_cal5.py [--tag k0882]
        [--results ansys_usermat/apdl/results/2026-10-cal5] [--out DIR]

Reads w8_<C>_cal_<tag>.json (stress, alpha of every element) and
w8_<C>_cal_<tag>_pm.json (last point-model state of the traced elements) as
written by run_calibrated.ps1 -Growth, for every condition that is present.
For the seed elements (seed_BIOFILM1) it gives the von Mises stress [Pa]
(deck units MPa, times 1e6), alpha - 1 (the UMAT variable), and the mean
volume fractions phi_i and living fractions phi_i psi_i. Writes
cal5_compare_<tag>.csv and .md and, unless --no-fig, cal5_compare_<tag>.png.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

CONDITIONS = ["CH", "CS", "DS", "DH"]
SPECIES = ["So", "An", "Vd", "Fn", "Pg"]


def seed_summary(run_json, pm_json):
    d = json.load(open(run_json))
    p = json.load(open(pm_json))
    seed = set(int(e) for e in d["seed_BIOFILM1"])
    a = d["all_stress"]
    elem = np.array(a["elem"]).astype(int)
    m = np.isin(elem, list(seed))
    seqv = np.array(a["seqv"], dtype=float)[m] * 1e6
    alpha = np.array(a["alpha"], dtype=float)[m]
    rows = [r for r in p["last_rows"] if int(r[0]) in seed and r[1] == "1"]
    if not rows:
        raise ValueError(f"{pm_json}: no point-model rows of the seed elements")
    c = p["columns"]
    i0 = c.index("gnew_1")
    g = np.array([[float(x) for x in r[i0:i0 + 12]] for r in rows])
    phi, psi = g[:, 0:5], g[:, 6:11]
    return {
        "n_seed": int(m.sum()),
        "seqv_mean": seqv.mean(), "seqv_max": seqv.max(),
        "alpha_mean": alpha.mean(),
        "phi": phi.mean(axis=0), "phipsi": (phi * psi).mean(axis=0),
        "theta_source": p.get("theta_source", ""),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="k0882")
    ap.add_argument("--results", default=str(Path(__file__).parent / "results" / "2026-10-cal5"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()
    res = Path(a.results)
    out = Path(a.out) if a.out else res
    sfx = f"_{a.tag}" if a.tag else ""

    found = {}
    for cnd in CONDITIONS:
        run = res / f"w8_{cnd}_cal{sfx}.json"
        pm = res / f"w8_{cnd}_cal{sfx}_pm.json"
        if run.exists() and pm.exists():
            found[cnd] = seed_summary(run, pm)
        else:
            print(f"{cnd}: not there yet ({run.name})")
    if not found:
        sys.exit("no runs found")

    ref = found.get("CH")
    head = ["condition", "n_seed", "seqv_mean_Pa", "seqv_max_Pa", "seqv_vs_CH",
            "alpha_minus_1_mean", "sum_phi", "sum_phipsi"] \
        + [f"phi_{s}" for s in SPECIES] + [f"phipsi_{s}" for s in SPECIES]
    table = []
    for cnd, s in found.items():
        table.append([cnd, s["n_seed"], s["seqv_mean"], s["seqv_max"],
                      s["seqv_mean"] / ref["seqv_mean"] if ref else float("nan"),
                      s["alpha_mean"], s["phi"].sum(), s["phipsi"].sum()]
                     + list(s["phi"]) + list(s["phipsi"]))

    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"cal5_compare{sfx}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        for r in table:
            w.writerow([r[0], r[1]] + [f"{x:.10e}" for x in r[2:]])

    md = [f"# Calibrated five-species runs, seed elements at T* = 1 (tag {a.tag or '-'})", "",
          "| condition | von Mises mean [Pa] | max [Pa] | / CH | alpha-1 mean | sum phi | sum phi psi | "
          + " | ".join(f"phi {s}" for s in SPECIES) + " |",
          "|---|" + "---|" * (6 + len(SPECIES))]
    for r in table:
        md.append(f"| {r[0]} | {r[2]:.4e} | {r[3]:.4e} | {r[4]:.3f} | {r[5]:.4e} | {r[6]:.4f} | {r[7]:.4f} | "
                  + " | ".join(f"{x:.4f}" for x in r[8:13]) + " |")
    md += ["", "Sources:"] + [f"- {c}: {s['theta_source']}" for c, s in found.items()]
    (out / f"cal5_compare{sfx}.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))

    if not a.no_fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import figstyle
        figstyle.apply()
        names = list(found)
        x = np.arange(len(names))
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
        ax1.bar(x, [found[c]["seqv_mean"] for c in names], color="0.4")
        ax1.set_xticks(x, names)
        ax1.set_ylabel("von Mises stress, seed mean [Pa]")
        bottom = np.zeros(len(names))
        for i, sp in enumerate(SPECIES):
            v = np.array([found[c]["phipsi"][i] for c in names])
            ax2.bar(x, v, bottom=bottom, label=sp)
            bottom += v
        ax2.set_xticks(x, names)
        ax2.set_ylabel(r"living fraction $\phi_i\psi_i$, seed mean")
        ax2.legend(frameon=False, fontsize=9, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.15))
        fig.tight_layout()
        fig.savefig(out / f"cal5_compare{sfx}.png", dpi=200)
        print(f"wrote {out / f'cal5_compare{sfx}.png'}")


if __name__ == "__main__":
    main()
