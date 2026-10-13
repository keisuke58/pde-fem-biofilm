# Heine et al. data: not in this repository

The measured data of Heine et al. (MHH/NIFE), for example the species-distribution
workbook `heine_species_distribution_biofilm.xlsx` and the raw qPCR and LIVE/DEAD
files, are available from the authors on request and are not redistributed here
(removed from the tree on 13 Oct 2026). Only aggregated summaries, at the level of
the published figures, are kept in git (`data/heine_volume/`).

Scripts that read the workbook (`plot_heine_phi_psi.py`, `plot_heine_composition.py`,
`ansys_usermat/apdl/clsm_phi_seed_reference.py`) expect a local copy at
`data/heine_species_distribution_biofilm.xlsx`; the tests that need it skip when it
is absent.
