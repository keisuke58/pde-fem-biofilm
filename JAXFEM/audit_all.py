"""
audit_all.py  —  IKM 修論品質オールインワン監査
================================================
5つのチェックを1コマンドで実行し、統合スコアを出す。

Run:
  python JAXFEM/audit_all.py           # 全チェック (RAG 含む)
  python JAXFEM/audit_all.py --quick   # RAG スキップ (高速, ≈3s)
  python JAXFEM/audit_all.py --eq      # equations のみ
  python JAXFEM/audit_all.py --fig     # figures のみ
  python JAXFEM/audit_all.py --reg     # regression のみ
  python JAXFEM/audit_all.py --thesis  # thesis_audit のみ
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
            print("     ℹ️  Known approximations acknowledged — thesis note at eq:alpha_proportional_split confirmed.")
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
    "eq":     ("Equations × Code",    run_equations),
    "fig":    ("FEM Figures",          run_figures),
    "reg":    ("FEM Regression",       run_regression),
    "thesis": ("Thesis LaTeX Audit",   run_thesis),
    "rag":    ("RAG Coverage",         run_rag),
}


def main():
    ap = argparse.ArgumentParser(
        description="IKM 修論品質オールインワン監査",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
          Examples:
            python JAXFEM/audit_all.py              # 全チェック
            python JAXFEM/audit_all.py --quick      # RAG スキップ (高速)
            python JAXFEM/audit_all.py --eq --fig   # 個別指定
        """))
    ap.add_argument("--quick",  action="store_true", help="Skip RAG (faster)")
    ap.add_argument("--eq",     action="store_true", help="Run equations check only")
    ap.add_argument("--fig",    action="store_true", help="Run figure comparison only")
    ap.add_argument("--reg",    action="store_true", help="Run regression guard only")
    ap.add_argument("--thesis", action="store_true", help="Run thesis audit only")
    ap.add_argument("--rag",    action="store_true", help="Run RAG coverage only")
    args = ap.parse_args()

    # Determine which sections to run
    flags = {k: getattr(args, k) for k in SECTIONS}
    any_flag = any(flags.values())

    if any_flag:
        to_run = {k: v for k, v in flags.items() if v}
    elif args.quick:
        to_run = {k: True for k in SECTIONS if k != "rag"}
    else:
        to_run = {k: True for k in SECTIONS}

    section_results = {}
    for key in SECTIONS:
        if not to_run.get(key):
            continue
        label, fn = SECTIONS[key]
        banner(label)
        out, counts = fn()
        # Print output (already captured but we want live output too)
        # Re-run live since redirect_stdout captured it
        section_results[key] = (out, counts)

    if len(section_results) > 1:
        rc = print_master_summary(section_results)
        sys.exit(rc)


if __name__ == "__main__":
    # Live output: don't redirect, just run each section directly
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick",  action="store_true")
    ap.add_argument("--eq",     action="store_true")
    ap.add_argument("--fig",    action="store_true")
    ap.add_argument("--reg",    action="store_true")
    ap.add_argument("--thesis", action="store_true")
    ap.add_argument("--rag",    action="store_true")
    args = ap.parse_args()

    flags = {k: getattr(args, k) for k in SECTIONS}
    any_flag = any(flags.values())

    if any_flag:
        to_run = {k: v for k, v in flags.items() if v}
    elif args.quick:
        to_run = {k: True for k in SECTIONS if k != "rag"}
    else:
        to_run = {k: True for k in SECTIONS}

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
