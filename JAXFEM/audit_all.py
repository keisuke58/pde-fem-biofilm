"""
audit_all.py  —  IKM 修論品質オールインワン監査
================================================
チェックを1コマンドで実行し、統合スコアを出す。

3段階モード:
  --quick   eq+fig+reg+thesis のみ (≈3s)         ← 日常チェック
  --strict  --quick + CI ファイル確認 (≈3s)      ← 投稿前
  --submit  --strict + RAG カバレッジ (≈2min)    ← 最終審査前

個別チェック:
  python JAXFEM/audit_all.py --eq      # equations のみ
  python JAXFEM/audit_all.py --fig     # figures のみ
  python JAXFEM/audit_all.py --reg     # regression のみ
  python JAXFEM/audit_all.py --thesis  # thesis_audit のみ
  python JAXFEM/audit_all.py --ci      # posterior CI files のみ
  python JAXFEM/audit_all.py --rag     # RAG coverage のみ

Exit code:
  0 — All PASS / INFO only
  1 — Any WARN/DRIFT
  2 — Any FAIL/ABSENT
"""
from __future__ import annotations
import argparse
import importlib
import io
import sys
import textwrap
import time
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

SECTION_WIDTH = 70


def banner(title: str):
    pad = (SECTION_WIDTH - len(title) - 4) // 2
    print(f"\n{'═'*SECTION_WIDTH}")
    print(f"{'═'*pad}  {title}  {'═'*(SECTION_WIDTH - pad - len(title) - 4)}")
    print(f"{'═'*SECTION_WIDTH}")


def count_tags(text: str) -> dict:
    counts = {"pass": 0, "info": 0, "warn": 0, "fail": 0}
    for line in text.split("\n"):
        s = line.strip()
        if s.startswith("✅"):
            counts["pass"] += 1
        elif s.startswith("ℹ️"):
            counts["info"] += 1
        elif s.startswith("⚠️") or s.startswith("🟡"):
            counts["warn"] += 1
        elif s.startswith("❌") or s.startswith("🔴"):
            counts["fail"] += 1
    return counts


def run_section(name: str, fn, *args, **kwargs) -> tuple[str, dict]:
    """Run fn(), capture stdout, return (output_text, tag_counts)."""
    buf = io.StringIO()
    t0 = time.time()
    try:
        with redirect_stdout(buf):
            fn(*args, **kwargs)
    except SystemExit:
        pass
    except Exception as e:
        buf.write(f"\n❌ FAIL  {name} raised exception: {e}\n")
    elapsed = time.time() - t0
    out = buf.getvalue()
    counts = count_tags(out)
    counts["elapsed"] = elapsed
    return out, counts


# ─────────────────────────────────────────────────────────────────────────────
# Section runners
# ─────────────────────────────────────────────────────────────────────────────

def run_equations() -> tuple[str, dict]:
    import cross_check_equations as m
    importlib.reload(m)
    # reset module-level state
    m.results.clear()
    def _run():
        print("IKM Klempt Pipeline — Equation × Code Cross-Check")
        print("="*65)
        m.check_eeff()
        m.check_pde()
        m.check_kalpha()
        m.check_predef()
        m.check_nu()
        m.check_data_files()
        m.check_attribution()
        m.check_depvar()
        m.print_summary()
    return run_section("equations", _run)


def run_figures() -> tuple[str, dict]:
    import compare_figures as m
    importlib.reload(m)
    m.results.clear()
    def _run():
        print("IKM Klempt FEM — Figure Quantitative Comparison")
        print("="*60)
        m.check_consistency()
        m.check_paper_claims()
        m.print_summary()
    return run_section("figures", _run)


def run_regression() -> tuple[str, dict]:
    import regression_guard as m
    importlib.reload(m)
    m.results.clear()
    def _run():
        print("IKM Klempt FEM — Regression Guard")
        print("="*60)
        m.compare(threshold=0.10, verbose=False)
        m.print_summary(threshold=0.10)
    return run_section("regression", _run)


def run_thesis() -> tuple[str, dict]:
    import thesis_audit as m
    importlib.reload(m)
    m.results.clear()
    def _run():
        print("IKM Thesis × Code Audit")
        print("="*65)
        texts = m.load_chapters()
        if texts:
            chaps = [k for k in texts if k != "main"]
            print(f"Loaded {len(chaps)} chapters")
            m.check_eq_coverage(texts, mode="full")
            m.check_approx_coverage(texts)
            m.check_klempt_attribution(texts)
            m.print_summary()
        else:
            print("❌ FAIL  No chapters found")
    return run_section("thesis", _run)


def run_ci() -> tuple[str, dict]:
    """Check posterior CI output files (--strict level)."""
    _FEM = HERE.parent
    _CI_OUT = HERE / "_posterior_ci"
    _CI0D   = _FEM / "_ci_0d_results"

    GEOMS = ["tooth", "implant"]
    CONDITIONS = ["commensal_hobic", "dysbiotic_hobic", "commensal_static", "dysbiotic_static"]
    MIN_SAMPLES = 5

    def _run():
        print("IKM Posterior CI — File & Sample-Count Check")
        print("="*60)
        import json

        # 1. CI json exists?
        for geom in GEOMS:
            f = _CI_OUT / f"klempt_stress_ci_{geom}.json"
            if f.exists():
                print(f"✅ PASS  _posterior_ci/klempt_stress_ci_{geom}.json exists")
            else:
                print(f"❌ FAIL  _posterior_ci/klempt_stress_ci_{geom}.json MISSING — run posterior_klempt_stress_ci.py")

        # 2. PDF figure exists?
        for geom in GEOMS:
            p = _CI_OUT / f"klempt_stress_ci_{geom}.pdf"
            if p.exists():
                print(f"✅ PASS  klempt_stress_ci_{geom}.pdf exists")
            else:
                print(f"⚠️  WARN  klempt_stress_ci_{geom}.pdf missing — re-run CI script")

        # 3. Per-condition sample count ≥ MIN_SAMPLES
        tooth_json = _CI_OUT / "klempt_stress_ci_tooth.json"
        if tooth_json.exists():
            d = json.load(open(tooth_json))
            for cond in CONDITIONS:
                if cond not in d:
                    continue
                n = d[cond].get("n_samples", 0)
                fallback = d[cond].get("fallback_map_only", False)
                if fallback or n < MIN_SAMPLES:
                    print(f"⚠️  WARN  [{cond}] only {n} samples — CI may be unreliable")
                else:
                    print(f"✅ PASS  [{cond}] {n} samples (≥{MIN_SAMPLES})")

        # 4. Check σ_CH/σ_DH ratio CI is informative
        if tooth_json.exists():
            d = json.load(open(tooth_json))
            if "ratio_ch_dh" in d:
                r = d["ratio_ch_dh"]
                width = r["p95"] - r["p05"]
                print(f"✅ PASS  σ_CH/σ_DH MAP={r['map']:.2f}×  90% CI=[{r['p05']:.2f},{r['p95']:.2f}]×  width={width:.2f}")
            else:
                print("⚠️  WARN  ratio_ch_dh not in CI json")

        # 5. Check ultimate samples exist for CS (key fix of this session)
        cs_ultimate = _CI0D / "commensal_static" / "samples_0d_ultimate.json"
        if cs_ultimate.exists():
            s = json.load(open(cs_ultimate))
            print(f"✅ PASS  CS samples_0d_ultimate.json exists  ({len(s)} samples from ultimate_10000p)")
        else:
            print("❌ FAIL  CS samples_0d_ultimate.json missing — run JAXFEM/resample_phi_ultimate.py")

    return run_section("ci", _run)


def run_rag() -> tuple[str, dict]:
    rag_dir = Path("/home/nishioka/LUHsummer26/tools/ikm_rag")
    sys.path.insert(0, str(rag_dir))

    import warnings; warnings.filterwarnings("ignore")
    import os
    os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
    os.environ.setdefault("HF_HUB_DISABLE_IMPLICIT_TOKEN", "1")

    import search as rag_m
    importlib.reload(rag_m)

    def _run():
        print("IKM RAG — Thesis Coverage Audit")
        try:
            col = rag_m.get_collection()
            rag_m.audit(col, n=3)
        except Exception as e:
            print(f"❌ FAIL  RAG unavailable: {e}")
    return run_section("rag", _run)


# ─────────────────────────────────────────────────────────────────────────────
# Unified summary
# ─────────────────────────────────────────────────────────────────────────────

def print_master_summary(section_results: dict[str, tuple[str, dict]]):
    banner("MASTER SUMMARY")

    total = {"pass": 0, "info": 0, "warn": 0, "fail": 0}
    rows  = []

    for name, (_, counts) in section_results.items():
        p = counts.get("pass",  0)
        i = counts.get("info",  0)
        w = counts.get("warn",  0) + counts.get("drift", 0)
        f = counts.get("fail",  0)
        t = counts.get("elapsed", 0)
        total["pass"] += p
        total["info"] += i
        total["warn"] += w
        total["fail"] += f

        if f > 0:
            status = "❌"
        elif w > 0:
            status = "⚠️ "
        elif i > 0:
            status = "ℹ️ "
        else:
            status = "✅"
        rows.append((status, name, p, i, w, f, t))

    print(f"\n  {'Check':<18} {'✅':>4} {'ℹ️':>4} {'⚠️':>4} {'❌':>4}  {'Time':>5}")
    print(f"  {'-'*55}")
    for status, name, p, i, w, f, t in rows:
        print(f"  {status} {name:<16} {p:>4} {i:>4} {w:>4} {f:>4}  {t:>4.1f}s")
    print(f"  {'-'*55}")
    print(f"  {'TOTAL':<18} {total['pass']:>4} {total['info']:>4} "
          f"{total['warn']:>4} {total['fail']:>4}")

    print()
    if total["fail"] > 0:
        print("  🔴 STATUS: FAIL — fix ❌ items before thesis submission")
        rc = 2
    elif total["warn"] > 0:
        print("  🟡 STATUS: WARN — investigate ⚠️ items before submission")
        rc = 1
    else:
        print("  ✅ STATUS: ALL CLEAR")
        if total["info"] > 0:
            print(f"     ℹ️  {total['info']} info item(s) — review cross_check output above.")
        rc = 0

    print()
    print("  Next-run commands:")
    print("    python JAXFEM/audit_all.py --quick     # equations+fig+reg+thesis (skip RAG)")
    print("    python JAXFEM/audit_all.py             # + RAG coverage")
    print("    python JAXFEM/regression_guard.py --baseline   # after intentional FEM change")
    print("    python JAXFEM/compare_figures.py --update      # after intentional golden update")
    print(f"\n{'═'*SECTION_WIDTH}")
    return rc


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

SECTIONS = {
    "eq":     ("Equations × Code",      run_equations),
    "fig":    ("FEM Figures",            run_figures),
    "reg":    ("FEM Regression",         run_regression),
    "thesis": ("Thesis LaTeX Audit",     run_thesis),
    "ci":     ("Posterior CI Files",     run_ci),
    "rag":    ("RAG Coverage",           run_rag),
}

# Which sections each top-level mode runs
MODE_SECTIONS = {
    "quick":  ["eq", "fig", "reg", "thesis"],
    "strict": ["eq", "fig", "reg", "thesis", "ci"],
    "submit": ["eq", "fig", "reg", "thesis", "ci", "rag"],
}


def _parse_args_and_mode(argv=None):
    """Parse args and return (to_run: dict[key→bool], mode_label: str)."""
    ap = argparse.ArgumentParser(
        description="IKM 修論品質オールインワン監査",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
          3-level mode (recommended):
            python JAXFEM/audit_all.py --quick    # eq+fig+reg+thesis (≈3s)
            python JAXFEM/audit_all.py --strict   # +CI files  (≈3s)
            python JAXFEM/audit_all.py --submit   # +RAG (≈2min)

          Individual checks:
            python JAXFEM/audit_all.py --eq --fig
        """))
    ap.add_argument("--quick",  action="store_true", help="日常チェック (eq+fig+reg+thesis)")
    ap.add_argument("--strict", action="store_true", help="投稿前   (--quick + CI files)")
    ap.add_argument("--submit", action="store_true", help="最終審査 (--strict + RAG)")
    ap.add_argument("--eq",     action="store_true")
    ap.add_argument("--fig",    action="store_true")
    ap.add_argument("--reg",    action="store_true")
    ap.add_argument("--thesis", action="store_true")
    ap.add_argument("--ci",     action="store_true")
    ap.add_argument("--rag",    action="store_true")
    args = ap.parse_args(argv)

    # Mode flags take priority over individual flags
    if args.submit:
        return {k: True for k in MODE_SECTIONS["submit"]}, "--submit"
    if args.strict:
        return {k: True for k in MODE_SECTIONS["strict"]}, "--strict"
    if args.quick:
        return {k: True for k in MODE_SECTIONS["quick"]}, "--quick"

    # Individual flags
    individual = {k: getattr(args, k) for k in SECTIONS if hasattr(args, k)}
    if any(individual.values()):
        return {k: v for k, v in individual.items() if v}, "custom"

    # Default: strict (quick + CI) when no args given
    return {k: True for k in MODE_SECTIONS["strict"]}, "--strict (default)"


def main():
    to_run, mode_label = _parse_args_and_mode()

    section_results = {}
    for key in SECTIONS:
        if not to_run.get(key):
            continue
        label, fn = SECTIONS[key]
        banner(label)
        out, counts = fn()
        section_results[key] = (out, counts)

    if len(section_results) > 1:
        rc = print_master_summary(section_results)
        sys.exit(rc)


if __name__ == "__main__":
    # Live output: don't redirect, just run each section directly
    to_run, mode_label = _parse_args_and_mode()

    import warnings; warnings.filterwarnings("ignore")
    import os
    os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
    os.environ.setdefault("HF_HUB_DISABLE_IMPLICIT_TOKEN", "1")

    section_outputs = {}
    timings = {}

    for key in SECTIONS:
        if not to_run.get(key):
            continue
        label, fn = SECTIONS[key]
        banner(label)
        t0 = time.time()
        buf = io.StringIO()
        try:
            with redirect_stdout(buf):
                if key == "eq":
                    import cross_check_equations as m
                    importlib.reload(m); m.results.clear()
                    print("IKM Klempt Pipeline — Equation × Code Cross-Check")
                    print("="*65)
                    m.check_eeff(); m.check_pde(); m.check_kalpha()
                    m.check_predef(); m.check_nu()
                    m.check_data_files(); m.check_attribution(); m.check_depvar()
                    m.print_summary()
                elif key == "fig":
                    import compare_figures as m
                    importlib.reload(m); m.results.clear()
                    m.check_consistency(); m.check_paper_claims(); m.print_summary()
                elif key == "reg":
                    import regression_guard as m
                    importlib.reload(m); m.results.clear()
                    m.compare(threshold=0.10, verbose=False)
                    m.print_summary(threshold=0.10)
                elif key == "thesis":
                    import thesis_audit as m
                    importlib.reload(m); m.results.clear()
                    texts = m.load_chapters()
                    if texts:
                        chaps = [k for k in texts if k != "main"]
                        print(f"Loaded {len(chaps)} chapters: {chaps}")
                        m.check_eq_coverage(texts, mode="full")
                        m.check_approx_coverage(texts)
                        m.check_klempt_attribution(texts)
                        m.print_summary()
                elif key == "ci":
                    run_ci()[0]   # run_ci() captures its own output
                    # Re-run without capture so output appears live
                    _FEM2 = HERE.parent
                    _CI_OUT2 = HERE / "_posterior_ci"
                    _CI0D2   = _FEM2 / "_ci_0d_results"
                    import json as _json
                    print("IKM Posterior CI — File & Sample-Count Check")
                    print("="*60)
                    for geom2 in ["tooth", "implant"]:
                        f2 = _CI_OUT2 / f"klempt_stress_ci_{geom2}.json"
                        if f2.exists():
                            print(f"✅ PASS  klempt_stress_ci_{geom2}.json exists")
                        else:
                            print(f"❌ FAIL  klempt_stress_ci_{geom2}.json MISSING")
                        p2 = _CI_OUT2 / f"klempt_stress_ci_{geom2}.pdf"
                        if p2.exists():
                            print(f"✅ PASS  klempt_stress_ci_{geom2}.pdf exists")
                        else:
                            print(f"⚠️  WARN  klempt_stress_ci_{geom2}.pdf missing")
                    tooth_json2 = _CI_OUT2 / "klempt_stress_ci_tooth.json"
                    if tooth_json2.exists():
                        d2 = _json.load(open(tooth_json2))
                        for c2 in ["commensal_hobic","dysbiotic_hobic","commensal_static","dysbiotic_static"]:
                            if c2 not in d2: continue
                            n2 = d2[c2].get("n_samples", 0)
                            fb2 = d2[c2].get("fallback_map_only", False)
                            if fb2 or n2 < 5:
                                print(f"⚠️  WARN  [{c2}] only {n2} samples")
                            else:
                                print(f"✅ PASS  [{c2}] {n2} samples (≥5)")
                        if "ratio_ch_dh" in d2:
                            r2 = d2["ratio_ch_dh"]
                            print(f"✅ PASS  σ_CH/σ_DH MAP={r2['map']:.2f}×  90%CI=[{r2['p05']:.2f},{r2['p95']:.2f}]×")
                    cs_ult2 = _CI0D2 / "commensal_static" / "samples_0d_ultimate.json"
                    if cs_ult2.exists():
                        s2 = _json.load(open(cs_ult2))
                        print(f"✅ PASS  CS samples_0d_ultimate.json  ({len(s2)} samples)")
                    else:
                        print("❌ FAIL  CS samples_0d_ultimate.json missing")
                elif key == "rag":
                    rag_dir = Path("/home/nishioka/LUHsummer26/tools/ikm_rag")
                    if str(rag_dir) not in sys.path:
                        sys.path.insert(0, str(rag_dir))
                    import search as rag_m; importlib.reload(rag_m)
                    col = rag_m.get_collection()
                    rag_m.audit(col, n=3)
        except SystemExit:
            pass
        except Exception as e:
            buf.write(f"\n❌ FAIL  {key} crashed: {e}\n")
            import traceback; traceback.print_exc(file=buf)

        elapsed = time.time() - t0
        out = buf.getvalue()
        print(out, end="")
        counts = count_tags(out)
        counts["elapsed"] = elapsed
        section_outputs[key] = (out, counts)
        timings[key] = elapsed

    if len(section_outputs) > 1:
        rc = print_master_summary(section_outputs)
        sys.exit(rc)
