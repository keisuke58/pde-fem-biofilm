# Experimental data

## `heine_species_distribution_biofilm.xlsx` — Heine in-vitro biofilm species distribution

CLSM/FISH-measured **relative species composition** of the 5-species oral biofilm
model (Heine et al., companion Nishioka–Heine paper), used to anchor the CLSM
composition φ in the pipeline (see the repo `README.md`).

- **Sheets** — one per state: `Commensal`, `Dysbiotic` (both *static, all cells*).
- **Species** (5 columns-blocks each):
  - Commensal: *S. oralis, A. naeslundii, V. dispar, F. nucleatum, P. gingivalis* (20709)
  - Dysbiotic: *S. oralis, A. naeslundii, V. parvula, F. nucleatum, P. gingivalis* (W83)
- **Rows** — `Tag` = timepoint in days: **1, 3, 6, 10, 15, 21**.
- **Values** — per-species relative abundance in **%** (each species block holds the
  multiple measurement columns / replicates; blanks = no measurement).

> Provenance: provided by the experimental side (Heine) as the raw
> species-distribution workbook. Kept verbatim; any derived/cleaned tables should
> be generated from this file, not edited in place.
