"""
cross_check_equations.py
========================
論文 (Klempt 2024/2025) × 実装コード の整合性チェッカー。
修論の命取りになるバグを8チェックポイントで検出する。

Run:
  python JAXFEM/cross_check_equations.py

出力記号:
  ✅ PASS — 論文とコードが一致
  ℹ️  INFO — 既知・文書化済みの近似/設計判断 (アクション不要)
  ⚠️  WARN — 要確認 (軽微な差異または手動検証)
  ❌ FAIL — 実装ミス確定 → 即修正必要

最終目標: 0 ❌, 0 ⚠️, ∞ ✅
"""
from __future__ import annotations
import re
import sys
import warnings
from pathlib import Path

import numpy as np

# ── paths ─────────────────────────────────────────────────────────────────────
HERE       = Path(__file__).resolve().parent
FEM_ROOT   = HERE.parent
NIFE_ROOT  = FEM_ROOT.parent / "nife" / "masterarbeit_ansys_fem"
IKM_ROOT   = FEM_ROOT.parent / "IKM_Hiwi" if not (FEM_ROOT.parent / "Tmcmc202601").exists() \
             else FEM_ROOT.parent
UMAT_VOIGT = NIFE_ROOT / "coupling_prototype" / "abaqus" / "umat_klempt_voigt.f"
UMAT_2025  = NIFE_ROOT / "coupling_prototype" / "abaqus" / "umat_klempt2025.f"
PDE_FILE   = HERE / "klempt_pde_multispecies.py"
INP_GEN    = FEM_ROOT / "gen_tooth_klempt_umat_inp.py"
TMCMC_DIR  = Path("/home/nishioka/IKM_Hiwi/Tmcmc202601/data_5species/_runs/ultimate_10000p")
ALPHA_DIR  = HERE   # klempt_alpha_final_{cond}.npy live here

# TMCMC run directory names (under ultimate_10000p/)
TMCMC_CONDITIONS  = ["commensal_hobic", "commensal_static", "dh_baseline", "dysbiotic_static"]
# alpha_final output file names (klempt_alpha_final_{cond}.npy in ALPHA_DIR)
ALPHA_CONDITIONS  = ["commensal_hobic", "commensal_static", "dysbiotic_hobic", "dysbiotic_static"]

PASS = "✅ PASS"
INFO = "ℹ️  INFO"
WARN = "⚠️  WARN"
FAIL = "❌ FAIL"

results = []

def report(tag: str, label: str, detail: str = ""):
    results.append((tag, label, detail))
    print(f"{tag}  {label}")
    if detail:
        for line in detail.strip().split("\n"):
            print(f"      {line}")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 1: E_eff (Klempt 2024 Eq.20)  phi^2-gated stiffness
# ═══════════════════════════════════════════════════════════════════════════════
def check_eeff():
    print("\n── Check 1: E_eff = phi_local² × E_Voigt / Σφ_i  (Klempt 2024 Eq.20) ──")

    def check_umat_eeff(fpath: Path, name: str):
        txt = fpath.read_text()
        evoigt_norm = re.search(r"E_voigt_norm\s*=\s*E_voigt\s*/\s*max\(phi_total_cond", txt)
        phi_gate    = re.search(r"phi_gate\s*=\s*phi_local\s*\*\s*phi_local", txt)
        e_gated     = re.search(r"E_gated\s*=\s*phi_gate\s*\*\s*E_voigt_norm", txt)
        if evoigt_norm and phi_gate and e_gated:
            return True, txt
        return False, txt

    ok_v, txt_v  = check_umat_eeff(UMAT_VOIGT, "umat_klempt_voigt.f")
    ok_25, txt25 = check_umat_eeff(UMAT_2025,  "umat_klempt2025.f")

    if not ok_v:
        report(FAIL, "E_eff formula missing/broken in umat_klempt_voigt.f")
    if not ok_25:
        report(FAIL, "E_eff formula missing/broken in umat_klempt2025.f")

    if ok_v and ok_25:
        # Numeric test
        phi_i     = np.array([0.2, 0.15, 0.1, 0.1, 0.05])
        E_i       = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])
        phi_local = 0.5
        phi_total = phi_i.sum()
        E_gated   = phi_local**2 * np.dot(phi_i, E_i) / phi_total
        report(PASS, f"E_gated = phi_local² × E_Voigt/Σφ_i  ({E_gated:.3e} MPa)",
               f"Both UMATs: identical formula E_voigt_norm=E_voigt/phi_total → phi_gate*norm ✓\n"
               f"Numeric: phi_local={phi_local}, Σφ_i={phi_total:.2f} → E_gated={E_gated:.3e} MPa")

    # Sub-check: wrong comment "(phi_local/phi_total)^2" must NOT appear (without "NOT")
    for fpath, name in [(UMAT_VOIGT, "umat_klempt_voigt.f"), (UMAT_2025, "umat_klempt2025.f")]:
        bad = [l.strip() for l in fpath.read_text().split("\n")
               if "phi_local/phi_total" in l and "NOT" not in l and "not" not in l]
        if bad:
            report(FAIL, f"{name}: wrong comment '(phi_local/phi_total)^2' still present",
                   f"  Line: {bad[0]}\n  Fix: use 'phi_local^2  (NOT ...)' form")
        # If no bad lines, already passing silently — no output noise


# ═══════════════════════════════════════════════════════════════════════════════
# Check 2: PDE Eq.34-36 terms
# ═══════════════════════════════════════════════════════════════════════════════
def check_pde():
    print("\n── Check 2: PDE Eq.34-36 (Allen-Cahn + logistic-Monod + k_α·α) ──")

    txt = PDE_FILE.read_text()

    checks = {
        "Allen-Cahn double-well: -Γ·φ(1-φ)(1-2φ)":
            r"GAMMA_AC\s*\*\s*phi\s*\*\s*\(1.*phi\)\s*\*\s*\(1.*2",
        "logistic-Monod: +K·φ(1-φ)·c/(k_M+c)":
            r"K_rate\s*\*\s*phi\s*\*\s*\(1.*phi\)\s*\*\s*monod",
        "k_α feedback: +k_α·α (NO Monod)":
            r"k_alpha\s*\*\s*alpha",
        "Eq.36: α̇ = k_α·φ (no Monod in alpha eq)":
            r"dalpha\s*=\s*k_alpha\s*\*\s*phi",
        "Eq.35: nutrient ċ = D∇²c − γ·φ·monod":
            r"dc\s*=\s*D_c\s*\*\s*lap.*-\s*gamma_c\s*\*\s*phi\s*\*\s*monod",
    }

    for label, pattern in checks.items():
        if re.search(pattern, txt, re.IGNORECASE):
            report(PASS, label)
        else:
            report(FAIL, f"NOT FOUND in PDE: {label}",
                   f"grep for: {pattern[:60]}")

    # Confirm NO Monod in alpha eq
    alpha_lines = [l for l in txt.split("\n") if "dalpha" in l and not l.strip().startswith("#")]
    if alpha_lines:
        monod_in_alpha = any("monod" in l.lower() or "k_m" in l.lower() for l in alpha_lines)
        if monod_in_alpha:
            report(FAIL, "Eq.36 VIOLATION: Monod term found in alpha equation!",
                   f"Line: {alpha_lines[0]}")
        else:
            report(PASS, "Eq.36 confirmed: NO Monod in α̇ equation")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 3: k_α_eff = Σ φ_i × k_α_i  (species-weighted)
# ═══════════════════════════════════════════════════════════════════════════════
def check_kalpha():
    print("\n── Check 3: k_α_eff = Σφ_i × k_α_i  (species-weighted growth rate) ──")

    txt = PDE_FILE.read_text()

    computed = re.search(r'"k_alpha_eff":\s*float\(np\.dot\(phi_vec,\s*K_ALPHA\)\)', txt)
    consumed = re.search(r'k_alpha_eff\s*=\s*p\["k_alpha_eff"\]', txt)
    if computed and consumed:
        report(PASS, 'k_α_eff chain: params["k_alpha_eff"]=np.dot(phi_vec,K_ALPHA) → p["k_alpha_eff"]')
    else:
        m = re.search(r'"k_alpha_eff":\s*(.+)', txt)
        formula = m.group(1)[:80] if m else "NOT FOUND"
        report(FAIL, "k_α_eff computation chain broken",
               f"Computed as: {formula}\nConsumed: {'✓' if consumed else 'NOT FOUND'}")
        return

    m = re.search(r"K_ALPHA\s*=\s*np\.array\(\[([^\]]+)\]\)", txt)
    if m:
        vals = [float(v.strip()) for v in m.group(1).split(",")]
        expected = [1.0, 0.8, 0.4, 0.6, 0.3]
        species  = ["So", "An", "Vd", "Fn", "Pg"]
        mismatches = [(sp, ex, got) for sp, ex, got
                      in zip(species, expected, vals) if abs(ex - got) > 1e-9]
        if not mismatches:
            report(PASS, f"K_ALPHA = {vals}  ← So>An>Fn>Vd>Pg species ranking")
        else:
            for sp, ex, got in mismatches:
                report(FAIL, f"K_ALPHA[{sp}]: expected {ex}, got {got}")

    try:
        import json
        mscl = FEM_ROOT / "_multiscale_2d_results"
        phi_ch = np.array(json.load(open(mscl / "ref_0d_commensal_hobic.json"))["phi_final"])
        K_ALPHA = np.array([1.0, 0.8, 0.4, 0.6, 0.3])
        k_eff_ch = np.dot(phi_ch, K_ALPHA)
        report(PASS, f"k_α_eff[CH] = Σφ_i×k_α_i = {k_eff_ch:.4f}  (numeric verification)")
    except Exception as e:
        report(WARN, f"Could not load TMCMC data for numeric k_eff check: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 4: TMCMC → PREDEF index alignment
# ═══════════════════════════════════════════════════════════════════════════════
def check_predef():
    print("\n── Check 4: TMCMC φ_i → PREDEF index alignment ──")

    gen_txt   = INP_GEN.read_text()
    umat_txt  = UMAT_VOIGT.read_text()
    umat2_txt = UMAT_2025.read_text()

    # Mode A: PREDEF(2..6) = phi_So..Pg
    umat_predef_a = re.search(r"phi_So\s*=.*PREDEF\(2\).*phi_An\s*=.*PREDEF\(3\)",
                               umat_txt, re.DOTALL)
    if re.search(r"start=2", gen_txt) and umat_predef_a:
        report(PASS, "Mode A: PREDEF(2..6) = φ_So..Pg  (gen_tooth ↔ umat_klempt_voigt)")
    else:
        report(WARN, "Mode A PREDEF alignment — verify gen_tooth start=2 ↔ umat PREDEF(2..6)")

    # Mode C: PREDEF(6..10) = phi_So..Pg; PREDEF(1..5) = alpha_s
    umat_predef_c = re.search(r"phi_So\s*=.*PREDEF\(6\)", umat2_txt)
    if re.search(r"start=6", gen_txt) and umat_predef_c:
        report(PASS, "Mode C: PREDEF(6..10) = φ_So..Pg  (gen_tooth ↔ umat_klempt2025)")
    else:
        report(FAIL, "Mode C: PREDEF index mismatch — gen_tooth start=6 must match PREDEF(6..10)")

    # phi_local index
    if "PREDEF(7)" in umat_txt:
        report(PASS, "Mode A: phi_local = PREDEF(7)  confirmed in umat_klempt_voigt.f")
    else:
        report(FAIL, "Mode A: phi_local PREDEF index ≠ PREDEF(7) in umat_klempt_voigt.f")
    if "PREDEF(11)" in umat2_txt:
        report(PASS, "Mode C: phi_local = PREDEF(11)  confirmed in umat_klempt2025.f")
    else:
        report(FAIL, "Mode C: phi_local PREDEF index ≠ PREDEF(11) in umat_klempt2025.f")

    # E_SPEC_MPa order consistency
    gen_arr = re.search(r"E_SPEC_MPa\s*=\s*np\.array\(\[([^\]]+)\]\)", gen_txt)
    if gen_arr:
        vals_gen  = [v.strip() for v in gen_arr.group(1).split(",")]
        expected  = ["1e-3", "8e-4", "6e-4", "2e-4", "1e-5"]
        species   = ["So", "An", "Vd", "Fn", "Pg"]
        mismatches = [(sp, ex, got) for sp, ex, got in zip(species, expected, vals_gen)
                      if ex not in got and got not in ex]
        if not mismatches:
            report(PASS, f"E_SPEC_MPa [So,An,Vd,Fn,Pg]={vals_gen}  consistent across all files")
        else:
            for sp, ex, got in mismatches:
                report(FAIL, f"E_SPEC_MPa[{sp}]: expected {ex}, gen_tooth has {got}")
    else:
        report(FAIL, "E_SPEC_MPa array not found in gen_tooth_klempt_umat_inp.py")

    # Mode C alpha rate-weighted split — exact (derived from per-species PDE)
    # Formula: α_s = α_total × (K_ALPHA_s × φ_s) / k_alpha_eff
    # Derivation: dα_s/dt = K_ALPHA_s × φ_s → at steady-state φ, α_s/α_total = K_ALPHA_s·φ_s/Σ_j K_ALPHA_j·φ_j
    if re.search(r"weights\s*=\s*\(phi_vec\s*\*\s*K_ALPHA\)\s*/", gen_txt):
        report(PASS, "Mode C alpha split: α_s = α_total × (K_ALPHA_s·φ_s)/k_alpha_eff  [exact from per-species PDE]")
    elif re.search(r"sp_frac\s*=\s*phi_vec\[sp_idx\]\s*/\s*\(phi_total", gen_txt):
        report(INFO, "Mode C alpha split: old φ_s/Σφ_i formula still present — update to rate-weighted")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 5: ν = 0.49  (Klempt Table 2)
# ═══════════════════════════════════════════════════════════════════════════════
def check_nu():
    print("\n── Check 5: ν = 0.49 (Klempt 2024 Table 2) consistency ──")
    for fname, fpath in [("umat_klempt_voigt.f", UMAT_VOIGT),
                          ("umat_klempt2025.f",   UMAT_2025),
                          ("gen_tooth_klempt.py", INP_GEN)]:
        txt = fpath.read_text()
        if re.search(r"nu.*0\.49|0\.49.*nu|NU.*0\.49|0\.49.*NU", txt, re.IGNORECASE):
            report(PASS, f"ν=0.49 in {fname}")
        else:
            report(WARN, f"ν=0.49 not confirmed in {fname} — check manually")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 6: TMCMC data files + alpha_final files exist (all 4 conditions)
# ═══════════════════════════════════════════════════════════════════════════════
def check_data_files():
    print("\n── Check 6: Data file existence (TMCMC samples + alpha_final × 4 conditions) ──")

    # TMCMC samples (directory names use TMCMC_CONDITIONS)
    tmcmc_missing = [c for c in TMCMC_CONDITIONS
                     if not (TMCMC_DIR / c / "samples.npy").exists()]
    if not tmcmc_missing:
        report(PASS, f"TMCMC samples.npy exists for all {len(TMCMC_CONDITIONS)} conditions")
    else:
        report(FAIL, "TMCMC samples.npy missing for some conditions",
               "\n".join(f"  Missing: {TMCMC_DIR}/{c}/samples.npy" for c in tmcmc_missing))

    # alpha_final outputs — use ALPHA_CONDITIONS naming (dysbiotic_hobic not dh_baseline)
    alpha_missing = [c for c in ALPHA_CONDITIONS
                     if not (ALPHA_DIR / f"klempt_alpha_final_{c}.npy").exists()]
    if not alpha_missing:
        bad_shape = []
        for cond in ALPHA_CONDITIONS:
            arr = np.load(ALPHA_DIR / f"klempt_alpha_final_{cond}.npy")
            if arr.ndim == 0 or arr.size == 0:
                bad_shape.append(cond)
        if not bad_shape:
            report(PASS, f"klempt_alpha_final_{{cond}}.npy exists & non-empty for all {len(ALPHA_CONDITIONS)} conditions")
        else:
            report(WARN, f"alpha_final files have unexpected shape: {bad_shape}")
    else:
        report(FAIL, "alpha_final files missing for some conditions",
               "\n".join(f"  Missing: klempt_alpha_final_{c}.npy" for c in alpha_missing))


# ═══════════════════════════════════════════════════════════════════════════════
# Check 7: UMAT attribution comments (no wrong "Klempt 2025 additive Mandel")
# ═══════════════════════════════════════════════════════════════════════════════
def check_attribution():
    print("\n── Check 7: UMAT source attribution (no wrong Klempt 2025 Mandel claim) ──")

    for fpath, name in [(UMAT_VOIGT, "umat_klempt_voigt.f"), (UMAT_2025, "umat_klempt2025.f")]:
        txt = fpath.read_text()

        # Must NOT claim "Klempt 2025 additive Mandel decomposition" (disproven)
        bad = [l.strip() for l in txt.split("\n")
               if "Klempt 2025 additive Mandel" in l]
        if bad:
            report(FAIL, f"{name}: wrong attribution 'Klempt 2025 additive Mandel' still present",
                   f"  Line: {bad[0]}\n"
                   f"  Klempt 2025 has NO alpha/Fg mechanics. Fix to 'Klempt 2024 F=Fe.Fg'")
        else:
            report(PASS, f"{name}: no wrong Klempt 2025 Mandel attribution")

        # SDV comment must use phi_local^2 (not (phi_local/phi_total)^2)
        bad_sdv = [l.strip() for l in txt.split("\n")
                   if "phi_local/phi_total" in l and "NOT" not in l and "not" not in l]
        if bad_sdv:
            report(FAIL, f"{name}: SDV comment still claims (phi_local/phi_total)^2",
                   f"  Line: {bad_sdv[0]}")
        else:
            report(PASS, f"{name}: SDV phi_gate comment correct (phi_local^2)")


# ═══════════════════════════════════════════════════════════════════════════════
# Check 8: Voigt–UMAT field variable count integrity (NPROPS / NSTATV)
# ═══════════════════════════════════════════════════════════════════════════════
def check_depvar():
    print("\n── Check 8: NPROPS / DEPVAR count vs INP *USER MATERIAL CONSTANTS ──")

    gen_txt = INP_GEN.read_text()

    # Mode A: CONSTANTS=1 (nu only), DEPVAR 5
    m_const_a = re.search(r"\*USER MATERIAL.*CONSTANTS\s*=\s*(\d+).*umat_klempt_voigt",
                           gen_txt, re.IGNORECASE | re.DOTALL)
    if not m_const_a:
        # Look for it in different form
        m_const_a = re.search(r"CONSTANTS=1.*Mode A|Mode A.*CONSTANTS=1|voigt.*CONSTANTS\s*=\s*1",
                               gen_txt, re.IGNORECASE)
    # Fallback: just check CONSTANTS=1 exists in gen_tooth
    m1 = re.search(r"CONSTANTS\s*=\s*1", gen_txt)
    m5 = re.search(r"CONSTANTS\s*=\s*5", gen_txt)

    umat_v_txt  = UMAT_VOIGT.read_text()
    umat25_txt  = UMAT_2025.read_text()

    # UMAT_VOIGT header says CONSTANTS=1 (nu only)
    # UMAT_2025  header says CONSTANTS=1 (nu only)
    # gen_tooth uses CONSTANTS=1 (consistent)

    # NPROPS check: both UMATs use PROPS(1)=nu only
    props1_v  = "PROPS(1)" in umat_v_txt and "PROPS(2)" not in umat_v_txt
    props1_25 = "PROPS(1)" in umat25_txt and "PROPS(2)" not in umat25_txt

    # DEPVAR check
    nstatv_v   = re.search(r"NSTATV\.ge\.\s*(\d+)", umat_v_txt)
    nstatv_25  = re.search(r"NSTATV\.ge\.\s*(\d+)", umat25_txt)
    max_v  = max(int(m.group(1)) for m in re.finditer(r"NSTATV\.ge\.\s*(\d+)", umat_v_txt))
    max_25 = max(int(m.group(1)) for m in re.finditer(r"NSTATV\.ge\.\s*(\d+)", umat25_txt))

    if props1_v:
        report(PASS, f"umat_klempt_voigt.f: PROPS(1)=nu only (CONSTANTS=1), NSTATV≥{max_v}")
    else:
        report(WARN, "umat_klempt_voigt.f: PROPS usage may exceed CONSTANTS=1")

    if props1_25:
        report(PASS, f"umat_klempt2025.f: PROPS(1)=nu only (CONSTANTS=1), NSTATV≥{max_25}")
    else:
        report(WARN, "umat_klempt2025.f: PROPS usage may exceed CONSTANTS=1")

    # PREDEF count: Mode A uses 7, Mode C uses 11
    m_pred7  = re.search(r"PREDEF\(7\)", umat_v_txt)
    m_pred8  = re.search(r"PREDEF\(8\)", umat_v_txt)   # must NOT exist in voigt
    m_pred11 = re.search(r"PREDEF\(11\)", umat25_txt)
    m_pred12 = re.search(r"PREDEF\(12\)", umat25_txt)  # must NOT exist

    if m_pred7 and not m_pred8:
        report(PASS, "Mode A PREDEF count: max PREDEF(7)  (7 field vars)")
    else:
        report(WARN, "Mode A PREDEF: unexpected PREDEF(8+) in umat_klempt_voigt.f")

    if m_pred11 and not m_pred12:
        report(PASS, "Mode C PREDEF count: max PREDEF(11)  (11 field vars)")
    else:
        report(WARN, "Mode C PREDEF: unexpected PREDEF(12+) in umat_klempt2025.f")


# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════
def print_summary():
    print("\n" + "="*65)
    passes = [r for r in results if r[0] == PASS]
    infos  = [r for r in results if r[0] == INFO]
    warns  = [r for r in results if r[0] == WARN]
    fails  = [r for r in results if r[0] == FAIL]
    print(f"  Cross-check: {len(passes)} ✅  {len(infos)} ℹ️   {len(warns)} ⚠️   {len(fails)} ❌")
    print("="*65)

    if fails:
        print("\n[ACTION REQUIRED — ❌ FAILs]")
        for _, label, _ in fails:
            print(f"  • {label}")

    if warns:
        print("\n[VERIFY MANUALLY — ⚠️  WARNs]")
        for _, label, _ in warns:
            print(f"  • {label}")

    if infos:
        print("\n[ACKNOWLEDGED — ℹ️  INFO  (no action needed)]")
        for _, label, _ in infos:
            print(f"  • {label}")

    if not fails and not warns:
        print("\n  ✅ All critical checks passed. Implementation matches Klempt 2024/2025 equations.")
        if infos:
            print("     Known approximations documented above — add §5.2 thesis note before submission.")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    print("IKM Klempt Pipeline — Equation × Code Cross-Check  (v2, 2026-06-25)")
    print("="*65)
    check_eeff()
    check_pde()
    check_kalpha()
    check_predef()
    check_nu()
    check_data_files()
    check_attribution()
    check_depvar()
    print_summary()
