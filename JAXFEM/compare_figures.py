"""
compare_figures.py
==================
FEM 出力の定量比較チェッカー。
tooth_klempt_comparison.json (ゴールデン) vs 現在の klempt_extract_*.csv を比較し、
±10 % 以内かどうかを PASS/FAIL で判定する。

Run:
  python JAXFEM/compare_figures.py               # 差分チェック
  python JAXFEM/compare_figures.py --update      # 現在値でゴールデンを更新
  python JAXFEM/compare_figures.py --ratio       # σ_CH/σ_DH 比のみ表示

出力記号:
  ✅ PASS  — 論文記載値と一致 (±10 %)
  ⚠️  DRIFT — 許容範囲外 → 原因調査が必要
  ❌ FAIL  — ファイルなし、または計算不可
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

# ── paths ─────────────────────────────────────────────────────────────────────
FEM_DIR    = Path(__file__).resolve().parent.parent
GOLDEN_F   = FEM_DIR / "tooth_klempt_comparison.json"
EXTRACT_DIR = FEM_DIR  # klempt_extract_{type}_{cond}.csv live here

CONDITIONS = ["commensal_hobic", "commensal_static", "dysbiotic_hobic", "dysbiotic_static"]
TYPES      = ["tooth", "implant"]
THRESHOLD  = 0.10   # 10 % tolerance

# Key paper claims (from HOLES_KLEMPT_MULTISPECIES.md, verified 2026-06-25)
#   tooth:   σ_CH / σ_DH (hobic) ≈ 6.1×  (CH has higher alpha → more growth → higher stress)
#   implant: σ_CH / σ_DH (hobic) ≈ 4.8×
PAPER_CLAIMS = {
    "tooth_CH_DH_ratio":   6.1,
    "implant_CH_DH_ratio": 4.8,
    "ratio_tol": 0.15,   # tighter: ±15 % on ratio
}

results = []

def report(tag: str, label: str, detail: str = ""):
    results.append((tag, label))
    sym = {"✅": "✅ PASS ", "⚠️": "⚠️  DRIFT", "❌": "❌ FAIL "}
    marker = next((v for k, v in sym.items() if tag.startswith(k)), tag)
    print(f"{marker}  {label}")
    if detail:
        for line in detail.strip().split("\n"):
            print(f"       {line}")


def load_extract(struct: str, cond: str) -> dict | None:
    f = EXTRACT_DIR / f"klempt_extract_{struct}_{cond}.csv"
    if not f.exists():
        return None
    import csv
    rows = list(csv.DictReader(open(f)))
    if not rows:
        return None
    col = "sigma_mises_MPa"
    if col not in rows[0]:
        return None
    vals = np.array([float(r[col]) for r in rows])
    alpha_col = "alpha" if "alpha" in rows[0] else None
    alpha_vals = np.array([float(r[alpha_col]) for r in rows]) if alpha_col else np.array([])
    return {
        "mises_max":  float(vals.max()),
        "mises_mean": float(vals.mean()),
        "alpha_max":  float(alpha_vals.max()) if len(alpha_vals) else None,
    }


def check_consistency():
    """Compare each condition × structure against golden JSON."""
    if not GOLDEN_F.exists():
        print(f"  Golden file not found: {GOLDEN_F}")
        print("  Run with --update to create it from current data.")
        return

    golden = json.load(open(GOLDEN_F))

    for struct in TYPES:
        print(f"\n── {struct.upper()} FEM: mises_max vs golden ──")
        for cond in CONDITIONS:
            curr = load_extract(struct, cond)
            if curr is None:
                report("❌", f"{struct}/{cond}: extract CSV missing")
                continue

            # Golden: tooth from JSON, implant computed at --update time
            gold_key = f"{struct}_{cond}_mises_max"
            gold_all = golden.get("_flat_golden", {})
            gold_v = gold_all.get(gold_key)

            if gold_v is None:
                # Try nested structure (tooth_klempt_comparison.json uses mode-nested form)
                for mode in ["A", "C"]:
                    v = golden.get(mode, {}).get(cond, {}).get("mises_max")
                    if v is not None:
                        gold_v = v
                        break

            if gold_v is None:
                report("⚠️", f"{struct}/{cond}: no golden reference — run --update",
                       f"current mises_max = {curr['mises_max']:.4e} MPa")
                continue

            ratio = curr["mises_max"] / gold_v
            delta_pct = abs(ratio - 1.0) * 100
            if delta_pct < THRESHOLD * 100:
                report("✅", f"{struct}/{cond}: mises_max = {curr['mises_max']:.4e} MPa  "
                             f"(Δ={delta_pct:.1f}% vs golden {gold_v:.4e})")
            else:
                report("⚠️", f"{struct}/{cond}: mises_max drifted {delta_pct:.1f}%",
                       f"  current  = {curr['mises_max']:.4e} MPa\n"
                       f"  golden   = {gold_v:.4e} MPa\n"
                       f"  ratio    = {ratio:.3f}  (allowed: 0.90–1.10)")


def check_paper_claims():
    """Verify σ_CH / σ_DH stress ratios against HOLES paper claims."""
    print("\n── σ_CH / σ_DH stress ratio (key paper claim) ──")

    for struct, claim_key, claim_val in [
        ("tooth",   "tooth_CH_DH_ratio",   PAPER_CLAIMS["tooth_CH_DH_ratio"]),
        ("implant", "implant_CH_DH_ratio",  PAPER_CLAIMS["implant_CH_DH_ratio"]),
    ]:
        ch = load_extract(struct, "commensal_hobic")
        dh = load_extract(struct, "dysbiotic_hobic")
        if ch is None or dh is None:
            report("❌", f"{struct}: extract CSVs missing for ratio check")
            continue
        if dh["mises_max"] < 1e-12:
            report("❌", f"{struct}: dysbiotic_hobic sigma is ~0 (check FEM run)")
            continue

        actual = ch["mises_max"] / dh["mises_max"]
        tol    = PAPER_CLAIMS["ratio_tol"]
        delta_pct = abs(actual - claim_val) / claim_val * 100

        if delta_pct < tol * 100:
            report("✅", f"{struct}: σ_CH/σ_DH = {actual:.2f}×  "
                         f"(paper ≈ {claim_val:.1f}×, Δ={delta_pct:.1f}%)")
        else:
            report("⚠️", f"{struct}: σ_CH/σ_DH = {actual:.2f}×  ≠  paper {claim_val:.1f}×",
                   f"  Δ = {delta_pct:.1f}%  (allowed ≤ {tol*100:.0f}%)\n"
                   f"  CH mises_max = {ch['mises_max']:.4e} MPa\n"
                   f"  DH mises_max = {dh['mises_max']:.4e} MPa\n"
                   f"  → Check UMAT E_SPEC_MPa, K_ALPHA, or alpha_final values")


def update_golden():
    """Overwrite golden values with current extract CSV data."""
    golden = {}
    if GOLDEN_F.exists():
        golden = json.load(open(GOLDEN_F))

    flat = {}
    for struct in TYPES:
        for cond in CONDITIONS:
            curr = load_extract(struct, cond)
            if curr is None:
                print(f"  SKIP {struct}/{cond}: CSV not found")
                continue
            flat[f"{struct}_{cond}_mises_max"]  = curr["mises_max"]
            flat[f"{struct}_{cond}_mises_mean"] = curr["mises_mean"]
            if curr["alpha_max"] is not None:
                flat[f"{struct}_{cond}_alpha_max"] = curr["alpha_max"]
            print(f"  Recorded {struct}/{cond}: mises_max={curr['mises_max']:.4e} MPa")

    golden["_flat_golden"] = flat
    with open(GOLDEN_F, "w") as f:
        json.dump(golden, f, indent=2)
    print(f"\n  Saved {len(flat)} values → {GOLDEN_F}")


def print_summary():
    passes  = sum(1 for r in results if r[0].startswith("✅"))
    drifts  = sum(1 for r in results if r[0].startswith("⚠️"))
    fails   = sum(1 for r in results if r[0].startswith("❌"))
    print(f"\n{'='*60}")
    print(f"  compare_figures: {passes} ✅  {drifts} ⚠️   {fails} ❌")
    print(f"{'='*60}")
    if drifts:
        print("\n[INVESTIGATE — ⚠️  DRIFTs]")
        for tag, label in results:
            if tag.startswith("⚠️"):
                print(f"  • {label}")
    if fails:
        print("\n[ACTION REQUIRED — ❌ FAILs]")
        for tag, label in results:
            if tag.startswith("❌"):
                print(f"  • {label}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="FEM figure quantitative comparison")
    ap.add_argument("--update", action="store_true",
                    help="Update golden_values.json from current extract CSVs")
    ap.add_argument("--ratio",  action="store_true",
                    help="Only check σ_CH/σ_DH paper claims")
    args = ap.parse_args()

    print("IKM Klempt FEM — Figure Quantitative Comparison")
    print("="*60)

    if args.update:
        print("\n[UPDATE MODE] Recording current CSV values as golden...")
        update_golden()
    elif args.ratio:
        check_paper_claims()
        print_summary()
    else:
        check_consistency()
        check_paper_claims()
        print_summary()
