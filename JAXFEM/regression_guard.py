"""
regression_guard.py
===================
FEM 結果のリグレッション防止ガード。
パラメータ変更後に FEM を実行した後これを走らせることで、
「静かに結果が変わった」事故を検出する。

Usage:
  python JAXFEM/regression_guard.py --baseline          # 現在値をゴールデンとして保存
  python JAXFEM/regression_guard.py                     # 現在値をゴールデンと比較
  python JAXFEM/regression_guard.py --threshold 0.05   # 閾値を 5 % に変更 (default 10 %)
  python JAXFEM/regression_guard.py --verbose           # 全メトリクスを表示

Golden file: FEM/regression_golden.json
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

FEM_DIR    = Path(__file__).resolve().parent.parent
GOLDEN_F   = FEM_DIR / "regression_golden.json"
EXTRACT_DIR = FEM_DIR

CONDITIONS = ["commensal_hobic", "commensal_static", "dysbiotic_hobic", "dysbiotic_static"]
TYPES      = ["tooth", "implant"]

# Metrics to track per condition × structure
METRICS = ["mises_max", "mises_mean", "mises_p95", "alpha_max", "alpha_mean", "E_gated_max"]

results = []

def report(tag: str, label: str, detail: str = ""):
    results.append((tag, label))
    print(f"{tag}  {label}")
    if detail:
        for line in detail.strip().split("\n"):
            print(f"      {line}")


def load_stats(struct: str, cond: str) -> dict | None:
    f = EXTRACT_DIR / f"klempt_extract_{struct}_{cond}.csv"
    if not f.exists():
        return None
    try:
        import csv
        rows = list(csv.DictReader(open(f)))
        if not rows:
            return None
        mises  = np.array([float(r["sigma_mises_MPa"]) for r in rows])
        alpha  = np.array([float(r["alpha"]) for r in rows]) if "alpha" in rows[0] else np.array([])
        egated = np.array([float(r["E_gated_MPa"]) for r in rows]) if "E_gated_MPa" in rows[0] else np.array([])
        stats = {
            "mises_max":   float(mises.max()),
            "mises_mean":  float(mises.mean()),
            "mises_p95":   float(np.percentile(mises, 95)),
        }
        if len(alpha):
            stats["alpha_max"]  = float(alpha.max())
            stats["alpha_mean"] = float(alpha.mean())
        if len(egated):
            stats["E_gated_max"] = float(egated.max())
        return stats
    except Exception as e:
        print(f"  ERROR reading {f}: {e}")
        return None


def save_baseline():
    golden = {"_version": "regression_guard_v1", "_metrics": METRICS, "data": {}}
    for struct in TYPES:
        golden["data"][struct] = {}
        for cond in CONDITIONS:
            stats = load_stats(struct, cond)
            if stats is None:
                print(f"  SKIP {struct}/{cond}: CSV missing")
                continue
            golden["data"][struct][cond] = stats
            print(f"  Recorded {struct}/{cond}:  "
                  f"mises_max={stats['mises_max']:.4e}  "
                  f"alpha_max={stats.get('alpha_max', 'N/A')}")

    with open(GOLDEN_F, "w") as f:
        json.dump(golden, f, indent=2)
    print(f"\n  Baseline saved → {GOLDEN_F}")
    n = sum(len(v) for v in golden["data"].values())
    print(f"  {n} condition records stored.")


def compare(threshold: float = 0.10, verbose: bool = False):
    if not GOLDEN_F.exists():
        print(f"  No baseline found at {GOLDEN_F}")
        print("  Run with --baseline to create it.")
        sys.exit(1)

    golden = json.load(open(GOLDEN_F))
    gold_data = golden.get("data", {})

    for struct in TYPES:
        print(f"\n── {struct.upper()} ──")
        gold_struct = gold_data.get(struct, {})

        for cond in CONDITIONS:
            curr = load_stats(struct, cond)
            gold = gold_struct.get(cond)

            if curr is None:
                report("❌", f"{struct}/{cond}: extract CSV missing")
                continue
            if gold is None:
                report("⚠️ ", f"{struct}/{cond}: no baseline entry — run --baseline")
                continue

            drifted = []
            lines   = []
            for metric in METRICS:
                g = gold.get(metric)
                c = curr.get(metric)
                if g is None or c is None:
                    continue
                if abs(g) < 1e-20:
                    continue
                pct = abs(c - g) / abs(g) * 100
                arrow = "↑" if c > g else "↓"
                line = (f"    {metric:<18}: {g:.4e} → {c:.4e}  "
                        f"{arrow}{pct:.1f}%{'  !!!' if pct > threshold*100 else ''}")
                lines.append(line)
                if pct > threshold * 100:
                    drifted.append((metric, g, c, pct))

            if not drifted:
                tag = "✅"
                label = (f"{struct}/{cond}: all metrics within ±{threshold*100:.0f}%")
                if verbose:
                    label += f"\n" + "\n".join(lines)
                report(tag, label)
            else:
                drift_list = ", ".join(f"{m}({p:.1f}%)" for m, _, _, p in drifted)
                report("⚠️ ", f"{struct}/{cond}: DRIFTED — {drift_list}",
                       "\n".join(lines))


def print_summary(threshold: float):
    passes = sum(1 for r in results if r[0] == "✅")
    drifts = sum(1 for r in results if r[0].startswith("⚠️"))
    fails  = sum(1 for r in results if r[0] == "❌")
    print(f"\n{'='*60}")
    print(f"  regression_guard (±{threshold*100:.0f}%): {passes} ✅  {drifts} ⚠️   {fails} ❌")
    print(f"{'='*60}")
    if drifts or fails:
        print("\n[INVESTIGATE]")
        for tag, label in results:
            if not tag.startswith("✅"):
                print(f"  {tag} {label.split(chr(10))[0]}")
        print(f"\n  Possible causes:")
        print(f"    1. UMAT parameters changed (E_SPEC_MPa, K_ALPHA, nu)")
        print(f"    2. alpha_final inputs changed (re-run klempt_pde_multispecies.py?)")
        print(f"    3. Mesh or boundary conditions modified")
        print(f"    4. Correct change → run --baseline to update golden")
    else:
        print("\n  All FEM metrics stable. No regression detected.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="FEM regression guard")
    ap.add_argument("--baseline",  action="store_true",
                    help="Save current results as regression baseline")
    ap.add_argument("--threshold", type=float, default=0.10,
                    help="Drift tolerance fraction (default 0.10 = 10%%)")
    ap.add_argument("--verbose",   action="store_true",
                    help="Show all metric values, not just drifted ones")
    args = ap.parse_args()

    print("IKM Klempt FEM — Regression Guard")
    print("="*60)

    if args.baseline:
        print("\n[BASELINE MODE] Saving current CSV stats as golden reference...")
        save_baseline()
    else:
        compare(threshold=args.threshold, verbose=args.verbose)
        print_summary(threshold=args.threshold)
