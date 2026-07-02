# Plan — Next Steps (2026-07-02)

Consolidated, prioritized roadmap for what comes after the 2026-06-29 report.
This ties together three existing threads that are currently scattered across
docs:

- the **research goals** ladder (Levels 1–4) — [research_goals_1_2.md](research_goals_1_2.md)
- the **rigor audit** limitations that gate publication — [VERIFICATION_SENSITIVITY_LIMITATIONS.md](VERIFICATION_SENSITIVITY_LIMITATIONS.md)
- the **remaining-tasks** table (A–F) at the end of `report_20260629.tex`

It supersedes the "next steps" scattered in the older subsystem plans
([FEM_PLAN.md](FEM_PLAN.md), [BENCHMARK_PLAN.md](BENCHMARK_PLAN.md),
[eigenstrain_theory_roadmap.md](eigenstrain_theory_roadmap.md)), which remain
valid as subsystem references but are dated (Feb 2026).

---

## 0. Where we are (snapshot)

| Component | Status |
|---|---|
| Klempt growth-stress core (`Fg=(1+α)I`, tangents, mesh convergence, dual-UMAT) | 🟢 verified (IKM standard) |
| Headline `σ_CH/σ_DH ≈ 6.4×` (early) with sensitivity bands | 🟢 reported with CI + limitations |
| Composition provenance (CLSM-anchored, TMCMC → A-matrix) | 🟢 corrected, documented |
| Extensions: AT2 phase-field, 2-chain visco UMAT, growth Jacobian, replicon gLV | 🟢 done 2026-06-29, **held for continuation** (not in thesis body) |
| Per-species stiffness `E_SPEC` | 🟡 assumed (needs species AFM) — reported as band |
| Per-condition growth `α` calibration | 🟡 magnitude-anchored only (noisy CLSM thickness) |
| Jaw-level / patient-specific pipeline (Level 2) | 🔴 not assembled as one system |

The **thesis core is done and defensible**. The open frontier is (a) closing the
two measurement-driven sensitivity gaps, and (b) turning the research scripts
into a jaw-level, uncertainty-aware pipeline — the "first paper" concept.

---

## 1. Three tracks and how to sequence them

There are three distinct destinations. They have different owners, timelines,
and risk. Do **not** interleave them on one branch.

1. **T1 — Thesis finalization (IKM submission).** Freeze the verified core; no
   new physics. Highest priority, shortest horizon.
2. **T2 — First paper (jaw-level mechano-ecological digital twin, Level 2).**
   The publishable next contribution; builds directly on the verified core.
3. **T3 — Muramatsu-lab continuation (research extensions A–F).** Longer-horizon
   research; the 2026-06-29 extensions are the seed.

Recommended order: **T1 → T2**, with **T3** items pulled forward only when they
unblock T2 (e.g. viscoelasticity if reviewers ask about time-dependence).

---

## 2. T1 — Thesis finalization (freeze the core)

Goal: a submittable, internally consistent thesis document with no loose claims.

- [ ] **Lock the headline framing.** Report the full `σ(t)` trajectory (Fig
      `sigma_trajectory_4cond`) as primary; quote `6.4×` only as the early-biofilm
      point with the S2/S3 bands attached. Cross-check every number against
      `VERIFICATION_SENSITIVITY_LIMITATIONS.md`.
- [~] **Verify all figures rebuild from committed code.** `JAXFEM/audit_all.py`
      gained an environment preflight: in an isolated checkout it now reports
      externally-sourced sections (Abaqus ODB extracts, sibling `../nife` UMATs,
      regression baseline, thesis/paper LaTeX + RAG) as SKIP instead of a
      misleading FAIL; `--quick` → ALL CLEAR (runnable subset), `--strict-env`
      forces the original behavior. A **true full ALL CLEAR still requires the
      author's complete workspace** (those inputs are prerequisites, not
      regressions) — run `audit_all.py --strict-env` there to confirm each thesis
      figure regenerates. See [PIPELINE.md](PIPELINE.md) § "Running the audit".
- [ ] **Confirm the DS-bug fix is reflected everywhere** the thesis cites DS
      (ref_0d, tooth/implant MAP, samples_0d) — see [DS_composition_fix.md](DS_composition_fix.md).
- [ ] **State scope explicitly.** In the limitations section, carry L1–L5
      verbatim in intent: `E_SPEC` contrast is assumed, `α` is magnitude-anchored,
      composition is CLSM-anchored (not TMCMC-inferred).
- [ ] **Decide extension inclusion.** Per the report, AT2 / visco / Jacobian /
      replicon are *held for continuation*. Keep them out of the thesis body;
      optionally reference as "future work" only.

**Acceptance:** thesis compiles, `audit_all.py` = ALL CLEAR, every quoted number
traces to a committed script, no claim exceeds what the audit supports.

---

## 3. T2 — First paper: jaw-level uncertainty-aware digital twin

This is the Level-2 target from `research_goals_1_2.md` (Fig1–Fig4 concept).
Working title (Option B): *"Towards a digital twin of periodontal biofilms:
TMCMC-calibrated mechano-ecological modeling on patient-specific jaws."*

The paper is a **method + jaw-level case study** — it does *not* claim clinical
outcome prediction. It shows how posterior parameter uncertainty propagates to
jaw-level stress/risk fields.

### T2.1 — Pipeline modularization (the enabling work) — ✅ implemented
The single biggest gap: the code runs as research scripts, not one pipeline.

- [x] Single entry point `pipeline.py`: `input (geometry + condition) → posterior
      → PDE α → stress CI → risk`, delegating each stage to the authoritative
      existing script (`run_posterior_pipeline.py`, `run_end_to_end_pipeline.py`,
      `posterior_klempt_stress_ci.py`, `risk_metric.py`) rather than rewriting.
- [x] Problem-specific config (`configs/*.json`) separated from generic stage
      machinery; stages SKIP gracefully when sibling-repo/Abaqus inputs are absent.
- [x] A config file per `{geom}×{condition}` (8 provided); switching is a
      parameter/CLI change, not a code edit.
- See [PIPELINE.md](PIPELINE.md). Runnable in an isolated checkout for the
  committed-input stages (`stress_ci` reuse + `risk`).

**Acceptance:** `python pipeline.py --config configs/tooth_commensal_hobic.json`
runs the committed-input stages through to a risk summary + figures.

### T2.2 — Patient-specific jaw geometry
- [ ] Establish a robust STL/CT → conformal mesh → BC path for a full jaw
      (Patient_1 dataset already referenced in `biofilm_conformal_tet.py` /
      `openjaw_*`). Make it re-runnable for ≥1 additional patient/tooth.
- [ ] Produce Fig1 (framework) and Fig3 (jaw-level stress, per-condition panels).

### T2.3 — Uncertainty propagation + risk metric — ✅ metric implemented
- [ ] Sample the TMCMC posterior → FEM ensemble → per-location mean + 95% credible
      band. Machinery exists: `posterior_uncertainty_propagation.py`,
      `run_posterior_abaqus_ensemble.py`, `aggregate_di_credible.py`,
      `JAXFEM/posterior_klempt_stress_ci.py`.
- [x] Clinically meaningful risk metric `P[σ > threshold]` implemented in
      `JAXFEM/risk_metric.py` (empirical posterior exceedance + 90% bootstrap CI +
      threshold-sweep survival curve). Unit-tested (`tests/test_risk_metric.py`,
      8 passing). See [PIPELINE.md](PIPELINE.md).
- [ ] Produce Fig4 (credible-band along a pocket line + jaw-surface risk map) —
      the per-location field version; the current metric is on the scalar
      peak-Mises posterior (per condition/geometry).

### T2.4 — Report generation
- [ ] Auto-generate a per-condition PDF/figure set (no full UI needed) —
      `generate_pipeline_summary.py` / `generate_paper_figures.py` are the seeds.

**Paper-readiness gate:** Fig1–Fig4 all regenerate from the modular pipeline for
≥1 patient and all 4 conditions, with uncertainty bands on the headline field.

---

## 4. Measurement gaps (gate the *quantitative* claims of T2)

These are the L1/L2 limitations. They cap how strong T2's numbers can be. None
are code tasks — flag them early to collaborators.

- [ ] **Species-specific AFM / rheology** → fixes the `E_SPEC` contrast (S2:
      3.7–12× band, ~half the ratio). Without it, T2 must keep the sensitivity
      band. Highest-value external measurement.
- [ ] **Clean per-condition thickness time series** → enables per-condition `α`
      calibration instead of magnitude-only (L2). CLSM z-profiles are currently
      noisy/non-monotonic and lack Static.

Until these land, T2 reports **ratios and bands**, not single absolute stresses —
consistent with the current framing.

---

## 5. T3 — Muramatsu-lab continuation (research extensions)

Seeded by the 2026-06-29 work. Priorities from the report's A–F table:

| ID | Item | Status | Next action |
|---|---|---|---|
| A | Growth Jacobian ∂σ/∂φᵢ | ✅ done | fold into a mechanistic sensitivity story |
| B | AT2 phase-field fracture | ✅ done | couple `t_crit` to the ecology (bridge exists) |
| D | 2-chain viscoelastic UMAT | ✅ done | run condition sweep if time-dependence matters |
| C′ | Replicon disordered-gLV | ✅ done | write up the small-community dysbiosis pattern |
| **C** | **State-dependent `A_ij` bifurcation** | 🔴 not started | **needs more data** — pursue after dataset grows |
| E | Neural PDE surrogate | 🔴 not started | low priority (speed only; direct solve works) |
| F | J-integral `G_c` post-processing | 🔴 not started | **B substitutes** — skip unless a reviewer asks |

**Most promising continuation thread:** the ecology↔mechanics bridge
(`replicon_at2_bridge.py`) — the finding that dysbiosis forces marginal ecological
stability yet slower mechanical fracture (ecological–mechanical trade-off) is a
genuinely novel, publishable mechanism. C (state-dependent interactions) is the
natural deepening but is data-gated.

---

## 6. Immediate next actions (this week)

Concrete, in priority order:

1. **T1 audit pass** — run `JAXFEM/audit_all.py`; fix any figure that doesn't
   regenerate from committed code. (Unblocks thesis freeze.)
2. **T2.1 pipeline consolidation** — inventory `run_end_to_end_pipeline.py` +
   `run_posterior_pipeline.py`, define the single config-driven entry point.
3. **T2.3 risk metric** — pick and implement one `P[σ > threshold]` metric on top
   of the existing posterior-propagation machinery.
4. **Measurement asks** — send the species-AFM and clean-thickness requests to
   collaborators now; they have the longest lead time.

---

## まとめ（日本語）

- **T1 修論**：検証済みコアを凍結。`σ(t)` 軌跡を主結果、`6.4×` は早期点＋感度帯で提示。
  拡張（AT2/粘弾性/Jacobian/replicon）は本体に入れず継続研究へ温存。
- **T2 第1報**：顎レベル・不確実性込みデジタルツイン（Level 2）。最大の作業は
  研究スクリプト群の**単一パイプライン化**（設定ファイル駆動）。事後分布 → FEM
  アンサンブル → リスク指標（P[σ>閾値]）→ Fig1–4 自動生成。
- **測定ギャップ**：菌種別 AFM（`E_SPEC` 対比）と清浄な厚さ時系列（`α` 校正）が
  定量主張の上限を決める。外部依頼はリードタイムが長いので今すぐ。
- **T3 継続研究**：生態↔力学ブリッジ（trade-off）が最有望。状態依存 `A_ij`（C）は
  データ増加後。E/F は低優先（F は B で代替可）。
</content>
</invoke>
