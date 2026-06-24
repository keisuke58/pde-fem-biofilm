"""
thesis_audit.py
===============
修論 LaTeX ソース → 実装コード の整合性ブリッジ。

やること:
  1. chapters/*.tex を全スキャン → 方程式ラベル / Klempt 引用 / 近似断り書きを抽出
  2. cross_check_equations.py の check 関数と対応付け
  3. 「論文で主張しているが実装検証がない」「実装で認めている近似が論文に書かれていない」を検出

Run:
  python JAXFEM/thesis_audit.py              # 全チェック
  python JAXFEM/thesis_audit.py --eqs        # 方程式ラベル一覧のみ
  python JAXFEM/thesis_audit.py --approx     # 近似チェックのみ (§5.2 確認用)

出力:
  ✅ OK      — 論文に記載あり + 実装で確認済み
  ⚠️  MISSING — 論文主張に実装検証なし (要 cross_check 追加)
  ℹ️  NOTE    — 実装近似が論文で明示されているか確認
  ❌ ABSENT  — 実装近似が論文に書かれていない (提出前に必須追記)
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

THESIS_DIR   = Path("/home/nishioka/LUHsummer26/1030_Masterarbeit")
CHAPTERS_DIR = THESIS_DIR / "chapters"
XCHECK       = Path(__file__).parent / "cross_check_equations.py"

results = []

def report(tag: str, label: str, detail: str = ""):
    results.append((tag, label))
    print(f"{tag}  {label}")
    if detail:
        for line in detail.strip().split("\n"):
            print(f"      {line}")


def load_chapters() -> dict[str, str]:
    texts = {}
    for tex in sorted(CHAPTERS_DIR.glob("*.tex")):
        texts[tex.stem] = tex.read_text()
    # Also load main.tex
    main = THESIS_DIR / "main.tex"
    if main.exists():
        texts["main"] = main.read_text()
    return texts


def extract_eq_labels(texts: dict[str, str]) -> dict[str, list[str]]:
    """Return {chapter: [label, ...]} for all \label{eq:...} in thesis."""
    out = {}
    for chap, text in texts.items():
        labels = re.findall(r"\\label\{(eq:[^}]+)\}", text)
        if labels:
            out[chap] = labels
    return out


def extract_eq_refs(texts: dict[str, str]) -> set[str]:
    """All \eqref{...} or \ref{eq:...} used in the thesis."""
    refs = set()
    for text in texts.values():
        refs.update(re.findall(r"\\eqref\{([^}]+)\}", text))
        refs.update(re.findall(r"\\ref\{(eq:[^}]+)\}", text))
    return refs


def extract_klempt_cites(texts: dict[str, str]) -> dict[str, list[str]]:
    """Where Klempt is cited (chapter → list of surrounding context)."""
    out = {}
    for chap, text in texts.items():
        matches = re.findall(r"[^\n]{0,60}\\cite[tp]?\{[^}]*[Kk]lempt[^}]*\}[^\n]{0,60}", text)
        if matches:
            out[chap] = matches
    return out


# ── Known equation mappings ──────────────────────────────────────────────────
# (thesis_label, cross_check_name, klempt_ref, description)
EQ_MAP = [
    ("eq:growth_split", "Check 1 / Check 4 / Check 7",
     "Klempt 2024 Eq.34-36",
     "F=Fe·Fg multiplicative growth split"),
    ("eq:neo_hookean", "Check 1",
     "Klempt 2024 Eq.21",
     "Neo-Hookean: σ = (lam·lnJe - mu)·I/Je + mu·be/Je"),
    ("eq:pde_strong", "Check 2",
     "Klempt 2024 Eq.34-36",
     "Allen-Cahn + logistic-Monod + k_α·α PDE"),
    ("eq:splitting", "Check 2",
     "Klempt 2024 Eq.35",
     "Reaction-diffusion splitting scheme"),
]

# ── Known approximations that MUST appear in thesis §5.x ────────────────────
REQUIRED_APPROX = [
    (
        "alpha_proportional_split",
        r"proportion|α_s\s*≈|alpha.*split|per.species.*growth.*approx",
        "§5.2",
        "Mode C per-species α_s = α_total × φ_s/Σφ_i  [Nishioka approx, no paper source]",
        "Must state: 'α_s ≈ α_total · φ_s/Σφ_i (proportional approximation)'",
    ),
    (
        "isotropic_growth_modelling_choice",
        r"isotropic.*growth.*modelling choice|volume.matched isotropy|lateral.*growth.*model",
        "§5.2",
        "F_g = J_g^{1/3}·I  isotropic growth for residual stress [explicit modelling choice]",
        "Already written at ch5 line 583-586 — verify present.",
    ),
    (
        "neo_hookean_fast_elastic",
        r"fast.elastic|short.time stress|visco.*relax|elastic.*idealisation",
        "§5.2",
        "Neo-Hookean = fast-elastic idealisation on seconds–minutes timescale",
        "Already written at ch5 line 600-603 — verify present.",
    ),
]


def check_eq_coverage(texts: dict, mode: str = "full"):
    print("\n── Equation label coverage ──")
    labels_by_chap = extract_eq_labels(texts)
    refs           = extract_eq_refs(texts)

    all_labels = {label for labels in labels_by_chap.values() for label in labels}

    if mode in ("eqs", "full"):
        print(f"\n  Defined equation labels ({len(all_labels)} total):")
        for chap, labs in labels_by_chap.items():
            for lab in labs:
                used = "✓ ref'd" if lab in refs else "  (unreferenced)"
                print(f"    \\label{{{lab}}}  [{chap}]  {used}")

    # Cross-check known equations
    print(f"\n  Known Klempt equations × cross_check mapping:")
    xcheck_txt = XCHECK.read_text() if XCHECK.exists() else ""
    for tex_label, cc_check, klempt_ref, desc in EQ_MAP:
        in_thesis = tex_label in all_labels
        in_xcheck = cc_check.split("/")[0].strip() in xcheck_txt if xcheck_txt else False
        if in_thesis and in_xcheck:
            report("✅", f"\\label{{{tex_label}}} present + {cc_check} verified",
                   f"  → {klempt_ref}: {desc}")
        elif in_thesis and not in_xcheck:
            report("⚠️ ", f"\\label{{{tex_label}}} in thesis but NOT in cross_check",
                   f"  → {klempt_ref}: {desc}\n  Add to cross_check_equations.py")
        else:
            report("ℹ️ ", f"\\label{{{tex_label}}} not yet defined in thesis",
                   f"  → Will need: {desc}")

    # Dangling references (eqref without label)
    dangling = refs - all_labels
    if dangling:
        print(f"\n  ⚠️  Dangling \\eqref{{}} (no matching \\label{{}}): {sorted(dangling)}")
    else:
        print(f"\n  ✅ No dangling \\eqref{{}} found")


def check_approx_coverage(texts: dict):
    """Verify that each known approximation is stated in the thesis."""
    print("\n── Approximation disclosure check ──")
    full_text = "\n".join(texts.values())

    for approx_id, pattern, section, description, guidance in REQUIRED_APPROX:
        found = bool(re.search(pattern, full_text, re.IGNORECASE))
        if found:
            report("✅", f"Approx '{approx_id}' stated in thesis",
                   f"  {description}")
        else:
            report("❌", f"Approx '{approx_id}' NOT found in thesis",
                   f"  Required at: {section}\n"
                   f"  Description: {description}\n"
                   f"  Guidance: {guidance}")


def check_klempt_attribution(texts: dict):
    """Scan for Klempt citations and confirm correct paper is credited."""
    print("\n── Klempt citation audit ──")
    klempt_cites = extract_klempt_cites(texts)

    if not klempt_cites:
        report("⚠️ ", "No Klempt citations found in any chapter",
               "  Klempt 2024 BMMB must be cited for F=Fe·Fg, Eq.20, Table 2")
        return

    total = sum(len(v) for v in klempt_cites.values())
    report("✅", f"Klempt cited {total} times across {len(klempt_cites)} chapters")

    # Warn if any context mentions "Mandel" + "2025" (the wrong attribution we fixed in code)
    for chap, ctxs in klempt_cites.items():
        for ctx in ctxs:
            if "2025" in ctx and ("Mandel" in ctx or "additive" in ctx.lower()):
                report("❌", f"Possible wrong attribution in {chap}:",
                       f"  Context: {ctx.strip()}\n"
                       f"  Klempt 2025 has NO Mandel/additive mechanics — use Klempt 2024")


def print_summary():
    ok      = sum(1 for r in results if r[0] == "✅")
    missing = sum(1 for r in results if r[0].startswith("⚠️"))
    absent  = sum(1 for r in results if r[0] == "❌")
    notes   = sum(1 for r in results if r[0] == "ℹ️")
    print(f"\n{'='*65}")
    print(f"  thesis_audit: {ok} ✅  {missing} ⚠️   {absent} ❌  {notes} ℹ️")
    print(f"{'='*65}")
    if absent:
        print("\n[MUST FIX BEFORE SUBMISSION — ❌ ABSENT]")
        for tag, label in results:
            if tag == "❌":
                print(f"  • {label}")
    if missing:
        print("\n[INVESTIGATE — ⚠️  MISSING]")
        for tag, label in results:
            if tag.startswith("⚠️"):
                print(f"  • {label}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--eqs",    action="store_true", help="Equation labels only")
    ap.add_argument("--approx", action="store_true", help="Approximation check only")
    ap.add_argument("--cites",  action="store_true", help="Klempt citation check only")
    args = ap.parse_args()

    print("IKM Thesis × Code Audit")
    print("="*65)

    texts = load_chapters()
    if not texts:
        print(f"ERROR: No .tex files found in {CHAPTERS_DIR}")
        sys.exit(1)

    found_chapters = [k for k in texts if k != "main"]
    print(f"Loaded {len(found_chapters)} chapters: {found_chapters}")

    if args.eqs:
        check_eq_coverage(texts, mode="eqs")
    elif args.approx:
        check_approx_coverage(texts)
    elif args.cites:
        check_klempt_attribution(texts)
    else:
        check_eq_coverage(texts, mode="full")
        check_approx_coverage(texts)
        check_klempt_attribution(texts)
        print_summary()
