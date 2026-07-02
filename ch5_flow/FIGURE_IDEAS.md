# Figure ideas — backlog & "already exists elsewhere" notes

Memo of candidate thesis/portfolio figures **not yet built as `ch5_flow/` TikZ**.
Key point (per review): most are *data/result* figures that **already have
generating scripts** in this repo or a sibling repo — reuse those rather than
redraw. Only a couple are genuinely new conceptual schematics.

Legend: ✅ built this session · ♻️ result figure, script/output likely already
exists (reuse) · 🆕 conceptual schematic worth drawing fresh in TikZ.

## Built (ch5_flow/ + umat_flow/)
✅ operator splitting, cross-diffusion FV, VE-vs-poro, 5-species Voigt UMAT,
Hamilton variational, data pipeline, CZM traction-separation, boundary
conditions, impl architecture, V&V convergence, time-scale separation,
multiphysics chain, **growth kinematics F=Fe·Fg**, **mixed-mode fracture (B-K)**;
plus the 4 `umat_flow/` UMAT algorithm flows.

## Backlog

| # | Idea | Status | Where it likely already lives / how to make it |
|---|---|---|---|
| 2 | σ_CH/σ_DH **sensitivity tornado** (growth floor ×3.4 × E_SPEC contrast → 3.7–12×) | ♻️ | `espec_sensitivity.py`, `plot_eeff_sensitivity.py`, `run_material_sensitivity_sweep.py`, `posterior_sensitivity_stress.py`, `plot_basin_sensitivity.py`; numbers in `VERIFICATION_SENSITIVITY_LIMITATIONS.md` (S2) + `espec_sensitivity_findings.md`. Tornado = small re-plot of these. |
| 3 | **TMCMC** tempering ladder + multimodal 15-D posterior (corner plot) | ♻️ | `generate_corner_plot_paper.py`, `plot_posterior_uncertainty.py`; posteriors in `../data_5species/_runs/*/samples.npy`, `../Tmcmc202601`. Corner plots almost certainly already produced for the TMCMC paper. |
| 4 | 5-species **gLV interaction network** (A-matrix), commensal vs dysbiotic | ♻️/🆕 | A-matrix stats from `JAXFEM/replicon_analysis.py`, `replicon_species_contrib.py`, `plot_species_competition.py` (posterior A in `../data_5species`). A *node-graph* rendering may be new — small TikZ/networkx on top of existing A. |
| 5 | **Uncertainty propagation** schematic + credible band / risk map | ♻️ | Result figures: `posterior_uncertainty_propagation.py`, `plot_posterior_uncertainty.py`, `aggregate_di_credible.py`, `run_posterior_abaqus_ensemble.py`; plus this session's `JAXFEM/risk_metric.py` / `risk_field.py`. A one-box *schematic* could complement, but the plots exist. |
| 7 | **V&V hierarchy pyramid** (code/solution verification, validation) | 🆕 | No existing asset — generic concept. Worth a fresh TikZ; complements `flow_vv_convergence`. |
| 8 | **Stress-relaxation / creep** master curves σ(t) | ♻️ | `generate_fig25_stress_relaxation.py`, `generate_fig26_creep.py`, `plot_visco_2ch.py`, `run_full_viscoelastic_analysis.py`. Already generated (fig25/fig26). |
| 9 | Anisotropy from ∇φ (transverse isotropy, DI lineage) | ♻️ | `run_aniso_comparison.py`, `fem_aniso_analysis.py`. |
| 11 | Percolation threshold / connectivity | ♻️ | `percolation_justification.tex` (already a LaTeX figure/derivation). |
| — | Growth Jacobian ∂σ/∂φᵢ | ♻️ | `JAXFEM/growth_jacobian.py`, `plot_growth_jacobian.py`. |
| — | 3-D tooth/implant **mesh + biofilm shell** geometry | ♻️ | Rendered (not TikZ): `fem_3d_visualize.py`, `plot_stress_3d*.py`, `odb_visualize.py`. |

## Suggested next actions
- If a figure is ♻️, **first look for the existing output/script** (this repo or
  `../data_5species` / `../Tmcmc202601` / `../nife`) before drawing anything.
- Genuinely new schematic worth adding: **7 (V&V pyramid)**, and optionally the
  **4 interaction-network graph** as a TikZ on top of the existing A-matrix.
