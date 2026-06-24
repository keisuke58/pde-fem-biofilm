"""
paper_audit.py
==============
投稿論文 2 本の LaTeX ソース監査。

対象:
  - hamilton_biofilm_nishioka.tex  (慶應課題研究, kadaikenkyu2026/)
  - nishioka_heine_paper.tex       (IKM/nife heine_paper/)

チェック項目:
  1. undefined citation  (.log から)
  2. undefined reference (.log から)
  3. missing figure      (.log から)
  4. TODO / placeholder  (tex 本文から)
  5. dangling \\eqref    (\\label と対応なし)
  6. 必須引用の存在確認  (論文別 required_cites)
  7. Klempt 帰属ミス検知 (Mandel+2025 の誤帰属)
  8. placeholder bib エントリ (bib 内の TODO: ...)

Run:
  python JAXFEM/paper_audit.py
"""
from __future__ import annotations
import re
from pathlib import Path

results = []

PASS = "✅"
WARN = "⚠️ "
FAIL = "❌"
INFO = "ℹ️ "


def report(tag: str, msg: str):
    results.append(tag)
    print(f"  {tag}  {msg}")


PAPERS = [
    {
        "name":    "hamilton_biofilm_nishioka",
        "tex":     Path("/home/nishioka/LUHsummer26/1050_Keio/kadaikenkyu2026/hamilton_biofilm_nishioka.tex"),
        "bib":     Path("/home/nishioka/LUHsummer26/1050_Keio/kadaikenkyu2026/reference.bib"),
        "log":     Path("/home/nishioka/LUHsummer26/1050_Keio/kadaikenkyu2026/hamilton_biofilm_nishioka.log"),
        "pdf":     Path("/home/nishioka/LUHsummer26/1050_Keio/kadaikenkyu2026/hamilton_biofilm_nishioka.pdf"),
        "required_cites": [
            "klempt2024",
            "klempt2025",
            "junker2021",
            "heine2025",
            "fritsch2025",
        ],
        "required_content": [
            (r"Hamilton|変分原理|Hamilton.*原理", "Hamilton 原理への言及"),
            (r"TMCMC|Transitional.*MCMC|posterior", "TMCMC/Bayes 推定への言及"),
            (r"biofilm|バイオフィルム", "バイオフィルムへの言及"),
        ],
    },
    {
        "name":    "klempt2024_hamilton_biofilm_JP",
        "tex":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_hamilton_biofilm_JP.tex"),
        "bib":     None,
        "log":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_hamilton_biofilm_JP.log"),
        "pdf":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_hamilton_biofilm_JP.pdf"),
        "required_cites": [],
        "required_content": [
            (r"Hamilton|変分原理", "Hamilton 原理への言及"),
            (r"biofilm|バイオフィルム|成長", "バイオフィルム成長への言及"),
            (r"Klempt|klempt", "Klempt への言及"),
        ],
    },
    {
        "name":    "klempt2024_slides",
        "tex":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_slides.tex"),
        "bib":     None,
        "log":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_slides.log"),
        "pdf":     Path("/home/nishioka/LUHsummer26/1030_Masterarbeit/notes/klempt2024_slides.pdf"),
        "required_cites": [],
        "required_content": [
            (r"Hamilton|変分原理", "Hamilton 原理への言及"),
            (r"biofilm|バイオフィルム|拡散", "バイオフィルム/拡散への言及"),
            (r"Klempt|klempt", "Klempt への言及"),
        ],
    },
    {
        "name":    "nishioka_heine_paper",
        "tex":     Path("/home/nishioka/IKM_Hiwi/nife/heine_paper/nishioka_heine_paper.tex"),
        "bib":     Path("/home/nishioka/IKM_Hiwi/nife/heine_paper/references_ikm.bib"),
        "log":     Path("/home/nishioka/IKM_Hiwi/nife/heine_paper/nishioka_heine_paper.log"),
        "pdf":     Path("/home/nishioka/IKM_Hiwi/nife/heine_paper/nishioka_heine_paper.pdf"),
        "required_cites": [
            "Heine2025PeriImplant",
            "Klempt2024ContinuumGrowth",
            "Klempt2025ContinuumBiofilm",
            "JunkerBalzani2021ExtendedHamilton",
            "Fritsch2025BayesianMicrofilms",
            "Dieckow2025InVivo",
            "ChingChen2007TMCMC",
            "Betz2016TMCMC",
        ],
        "required_content": [
            (r"Hamilton|variational.*principle|Lagrangian", "Hamilton principle"),
            (r"TMCMC|transitional.*MCMC|posterior", "TMCMC"),
            (r"dysbiosis|commensal|peri.implant", "dysbiosis / peri-implant"),
            (r"interaction.*matrix|A_\{ij\}|\\mathbf\{A\}", "interaction matrix A_ij"),
        ],
    },
]


def check_log(paper: dict):
    """Parse .log file for undefined cites, refs, missing figures."""
    log = paper["log"]
    name = paper["name"]
    if not log.exists():
        report(WARN, f"[{name}] .log not found — build first")
        return

    text = log.read_text(errors="replace")

    # Undefined citations
    undef_cites = sorted(set(re.findall(
        r"Citation `([^']+)' on page \d+ undefined", text)))
    if undef_cites:
        for c in undef_cites:
            report(FAIL, f"[{name}] undefined \\cite{{{c}}}")
    else:
        report(PASS, f"[{name}] no undefined citations")

    # Undefined labels/refs
    undef_refs = sorted(set(re.findall(
        r"Reference `([^']+)' on page \d+ undefined", text)))
    if undef_refs:
        for r in undef_refs:
            report(WARN, f"[{name}] undefined \\ref/{{{r}}}")
    else:
        report(PASS, f"[{name}] no undefined references")

    # Missing figures
    missing_figs = sorted(set(re.findall(
        r"File `([^']+)' not found", text)))
    if missing_figs:
        for f in missing_figs:
            report(INFO, f"[{name}] missing figure: {f}")
    else:
        report(PASS, f"[{name}] all included files found")

    # Overfull hboxes
    n_overfull = len(re.findall(r"Overfull \\hbox", text))
    if n_overfull > 10:
        report(WARN, f"[{name}] {n_overfull} Overfull \\hbox warnings")
    elif n_overfull > 0:
        report(INFO, f"[{name}] {n_overfull} Overfull \\hbox (minor)")
    else:
        report(PASS, f"[{name}] no Overfull \\hbox")


def check_tex(paper: dict):
    """Check tex source for TODOs, dangling eqrefs, required cites/content."""
    tex = paper["tex"]
    name = paper["name"]
    if not tex.exists():
        report(FAIL, f"[{name}] tex file not found: {tex}")
        return

    text = tex.read_text(errors="replace")

    # TODO / FIXME / placeholder
    todos = re.findall(r"(?i)(TODO|FIXME|PLACEHOLDER|XX+|MISSING|TBD|\\textbf\{TODO\})[^\n]{0,60}", text)
    if todos:
        for t in todos[:5]:
            report(WARN, f"[{name}] TODO/placeholder: {t[:60].strip()}")
        if len(todos) > 5:
            report(WARN, f"[{name}] … {len(todos)-5} more TODO/FIXME/PLACEHOLDER")
    else:
        report(PASS, f"[{name}] no TODO/FIXME/PLACEHOLDER markers")

    # Dangling \eqref: refs without matching \label
    labels = set(re.findall(r"\\label\{([^}]+)\}", text))
    eqrefs = set(re.findall(r"\\eqref\{([^}]+)\}", text))
    dangling = eqrefs - labels
    if dangling:
        for d in sorted(dangling):
            report(WARN, f"[{name}] dangling \\eqref{{{d}}} (no matching \\label)")
    else:
        report(PASS, f"[{name}] no dangling \\eqref")

    # Required citations present
    for cite_key in paper.get("required_cites", []):
        pattern = r"\\cite[tp]?\{[^}]*" + re.escape(cite_key) + r"[^}]*\}"
        if re.search(pattern, text):
            report(PASS, f"[{name}] \\cite{{{cite_key}}} present")
        else:
            report(FAIL, f"[{name}] \\cite{{{cite_key}}} NOT FOUND — required citation missing")

    # Required content patterns
    for pattern, desc in paper.get("required_content", []):
        if re.search(pattern, text, re.IGNORECASE):
            report(PASS, f"[{name}] content present: {desc}")
        else:
            report(WARN, f"[{name}] content possibly missing: {desc}")

    # Klempt 2025 wrong attribution (Mandel / additive)
    bad = re.findall(
        r"[^\n]{0,40}(?:Klempt.*2025|2025.*Klempt)[^\n]{0,40}(?:Mandel|additive decomp)[^\n]{0,40}",
        text, re.IGNORECASE)
    if bad:
        for b in bad:
            report(FAIL, f"[{name}] Klempt 2025 wrong attribution: {b.strip()[:80]}")
    else:
        report(PASS, f"[{name}] no Klempt 2025 Mandel mis-attribution")


def check_bib(paper: dict):
    """Check bib for placeholder TODO entries."""
    bib = paper["bib"]
    name = paper["name"]
    if bib is None:
        return
    if not bib.exists():
        report(WARN, f"[{name}] bib file not found: {bib}")
        return

    text = bib.read_text(errors="replace")
    todos = re.findall(r"TODO[^\n]{0,80}", text)
    if todos:
        for t in todos:
            report(INFO, f"[{name}] bib TODO: {t[:70].strip()}")
    else:
        report(PASS, f"[{name}] no TODO placeholders in bib")

    # Entries that exist in registry (dynamic_registry.json) but have no PDF are skipped at runtime
    # just count total entries
    n_entries = len(re.findall(r"^@", text, re.MULTILINE))
    report(INFO, f"[{name}] bib: {n_entries} entries")


def check_pdf(paper: dict):
    """Check PDF exists and has reasonable size."""
    pdf = paper["pdf"]
    name = paper["name"]
    if pdf.exists():
        kb = pdf.stat().st_size // 1024
        if kb < 50:
            report(WARN, f"[{name}] PDF very small ({kb} KB) — may be empty/corrupt")
        else:
            report(PASS, f"[{name}] PDF exists ({kb} KB)")
    else:
        report(FAIL, f"[{name}] PDF not found — build first")


def run():
    for paper in PAPERS:
        name = paper["name"]
        print(f"\n{'─'*65}")
        print(f"  Paper: {name}")
        print(f"{'─'*65}")
        check_pdf(paper)
        check_log(paper)
        check_tex(paper)
        check_bib(paper)

    print(f"\n{'─'*65}")
    n_pass = sum(1 for r in results if r == PASS)
    n_warn = sum(1 for r in results if r == WARN)
    n_info = sum(1 for r in results if r == INFO)
    n_fail = sum(1 for r in results if r == FAIL)
    print(f"  Paper Audit Summary: {n_pass} ✅  {n_info} ℹ️   {n_warn} ⚠️   {n_fail} ❌")


if __name__ == "__main__":
    run()
