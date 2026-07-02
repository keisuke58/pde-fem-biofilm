# Pipeline & risk metric (T2.1 + T2.3)

Implements two items from [PLAN_NEXT.md](PLAN_NEXT.md):

- **T2.1** — a single, config-driven entry point (`pipeline.py`) that unifies the
  previously separate `run_posterior_pipeline.py` and `run_end_to_end_pipeline.py`
  behind one orchestrator, separating problem-specific config (JSON under
  `configs/`) from the generic stage machinery.
- **T2.3** — a clinical risk metric `P[σ > threshold]` (`JAXFEM/risk_metric.py`)
  as the terminal stage, computed from the posterior stress samples.

---

## Quick start

```bash
# One condition, all stages (skips stages whose inputs aren't in this checkout):
python pipeline.py --config configs/tooth_commensal_hobic.json

# Just the risk metric on the committed posterior CI:
python pipeline.py --config configs/tooth_commensal_hobic.json --stages risk

# Standalone risk metric (both geometries):
python JAXFEM/risk_metric.py --geom tooth --threshold-kpa 5.0
python JAXFEM/risk_metric.py --geom implant

# See the plan without running anything:
python pipeline.py --config configs/tooth_commensal_hobic.json --dry-run
python pipeline.py --list-stages
```

Requires `numpy` (+ `matplotlib` for figures). Tests: `python -m pytest tests/test_risk_metric.py`.

---

## Stages

`pipeline.py` runs an ordered subset of four stages. Each **checks its
prerequisites first** and SKIPs (with a reason) rather than crashing when inputs
are absent — so the pipeline is runnable in an isolated checkout for the stages
whose inputs are committed.

| Stage | Delegates to | Produces | Inputs (and where) |
|---|---|---|---|
| `posterior` | `run_posterior_pipeline.py` | posterior FEM ensemble | TMCMC run dir (`../data_5species/_runs`, sibling repo) |
| `forward` | `run_end_to_end_pipeline.py` | PDE α → eigenstrain → INP | `theta_MAP.json` (sibling repo) |
| `stress_ci` | `JAXFEM/posterior_klempt_stress_ci.py` | `_posterior_ci/klempt_stress_ci_{geom}.json` | 0D samples (`_ci_0d_results`, gitignored) |
| `risk` | `JAXFEM/risk_metric.py` | `risk_summary_{geom}.json` + figures | the CI json above (**committed**) |

The `risk` stage runs fully in a fresh checkout because the posterior stress
samples are committed in `JAXFEM/_posterior_ci/`.

### `stress_ci` is reuse-by-default

The committed CI json already carries the full per-condition posterior samples.
Recomputing needs the upstream 0D samples for **all four** conditions; a partial
set would silently collapse the missing conditions to MAP-only and clobber the
committed json. So `stress_ci` reuses the committed json unless the config sets
`"stress_ci": {"recompute": true}` **and** every condition's 0D samples are
present.

---

## Config format (`configs/*.json`)

Problem-specific settings live entirely in the config; the code is generic.

```json
{
  "condition": "commensal_hobic",
  "geom": "tooth",
  "stages": ["posterior", "forward", "stress_ci", "risk"],
  "risk": {"threshold_kpa": 5.0, "reference_kpa": [2.5, 5.0, 10.0]},
  "posterior": {"n_samples": 20, "nx": 15, "ny": 15, "nz": 15},
  "forward": {"quick": false, "no_inp": true},
  "stress_ci": {"recompute": false}
}
```

Eight examples are provided (`{tooth,implant} × {CH,CS,DS,DH}`). CLI flags
`--condition`, `--geom`, `--stages` override the config. A `pipeline_manifest.json`
recording what ran/skipped is written to `_pipeline/{geom}_{condition}/`.

---

## Risk metric — `P[σ > threshold]`

`JAXFEM/risk_metric.py` turns the per-condition posterior stress distribution
into an exceedance probability: the fraction of the posterior with peak von
Mises stress above a detachment-relevant threshold `τ`.

- **Threshold is not load-bearing.** Because biofilm cohesive/detachment strength
  is not a well-established constant, the full survival curve `P[σ > τ]` is the
  primary read-out (`risk_survival_{geom}.png`); point values are reported at
  several reference thresholds (default 2.5 / 5 / 10 kPa).
- **Uncertainty.** Each exceedance probability carries a 90% bootstrap CI
  (per-condition posterior samples are small, n ≈ 22–51). Degenerate
  (near-constant) conditions collapse to a point interval — the honest answer.
- **Units.** Stresses are Abaqus MPa; thresholds are in kPa (1 kPa = 1e-3 MPa).

Example (tooth, τ = 5 kPa): CH (commensal-HOBIC, early So-dominant biofilm)
carries the highest growth-stress exceedance (P ≈ 1.0), consistent with the
headline `σ_CH/σ_DH` result; DS ≈ 0. Outputs go to `JAXFEM/_risk/` (gitignored;
regenerate on demand — the summary json is the durable artifact).

---

## Running the audit in an isolated checkout

`JAXFEM/audit_all.py` gained an **environment preflight**: sections whose inputs
live outside this repo (Abaqus ODB extracts, the sibling `../nife` UMAT sources,
a regression baseline, and the thesis/paper LaTeX + RAG index in the author's
workspace) are reported as **SKIP** (neutral) instead of a misleading ❌ FAIL.

```bash
python JAXFEM/audit_all.py --quick        # SKIPs external sections → ALL CLEAR (runnable subset)
python JAXFEM/audit_all.py --ci           # runs: posterior CI files are committed
python JAXFEM/audit_all.py --strict-env   # force every section (original fail-loud behavior)
```

A **true full ALL CLEAR** still requires the author's complete workspace: the
`../nife` UMAT sources (eq), Abaqus ODB extract CSVs (fig), a
`regression_golden.json` baseline (reg), the thesis chapters and paper LaTeX
(thesis/papers), and the gitignored 0D ultimate samples (one `ci` sub-check).
Those are prerequisites, not regressions.
