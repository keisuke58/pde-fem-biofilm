# Experimental data for comparison (search of 4 Nov 2026)

Collected for the Keio continuation (KEIO_PLAN.ja.md 0b, idea 4: compare the model with measurements). Values are as
reported in the abstracts or snippets; check the full texts before quoting them. The model stresses scale linearly with
the stiffness E (thesis: E = 10 Pa from Klempt et al. 2024, stresses of order 1e-4 Pa).

## Saved in this folder (CC BY 4.0; credit in THIRD_PARTY.md)
| file | paper | data taken from the full text |
|---|---|---|
| `Pattem2018_oral_biofilm_AFM_SciRep.pdf/.txt` | Pattem et al. 2018, Sci. Rep. 8:5691, doi:10.1038/s41598-018-23798-1 | microcosm oral biofilms on hydroxyapatite, AFM (Hertz, nu = 0.5). Young's modulus, two populations: **day 3**: 14.35 ± 1.75 to 41.05 ± 0.97 kPa (0.1 % sucrose, nutrient-poor), 0.55 ± 0.02 to 2.57 ± 0.17 kPa (5 % sucrose, nutrient-rich); **day 5**: 1.17 ± 0.08 to 46.07 ± 1.87 kPa and 0.56 ± 0.06 to 6.66 ± 0.46 kPa. No significant change between day 3 and 5 within a group. More sucrose: softer, more adhesive. |
| `Granicher2021_CHX_multispecies_FrontMicrobiol.pdf/.txt` | Gränicher et al. 2021, Front. Microbiol. 12:741863 | six-species supragingival biofilm (Zurich model: S. mutans OMZ 918, S. oralis OMZ 607, A. oris OMZ 745, V. dispar OMZ 493, F. nucleatum, C. albicans), 64 h; chlorhexidine 0.05, 0.1 and 0.2 %, with and without DNase I / proteinase K pretreatment; CFU per species (Fig. 1) and CLSM live/dead. A dose-response end point per species, not a time course. |
| `Geisel2022_biofilm_wrinkling_eLife.txt` (text only, PDF 12 MB: https://elifesciences.org/articles/76027) | Geisel et al. 2022, eLife 11:e76027 | P. aeruginosa in a microfluidic channel; growth under lateral confinement -> compressive stress -> buckling-delamination (critical stress after Hutchinson et al. 1992); wall shear stress 0.1 Pa, far below the biofilm shear modulus. |

Use for the model: Pattem gives the stiffness range (the model's 10 Pa is 50 to 4000 times lower; the stress scales with E);
the sucrose dependence is a natural test of a composition- or nutrient-dependent stiffness (idea 2 of KEIO_PLAN 0b).
Gränicher gives per-species survival under an antiseptic dose, a candidate end point for the point model's antibiotic term.

## Stiffness of oral biofilms (from the search)
| source | what | value | note |
|---|---|---|---|
| Pattem et al. 2018 (already cited in the thesis) | oral biofilm, AFM | 0.6 to 41 kPa | depends on the sucrose supply |
| Vinogradov et al. 2004, *Biofilms* | *S. mutans*, shear rheometer, rotating-disk reactor | creep, G', G'', viscosity | viscoelastic fluid-like behaviour |
| Winston et al. (Montana State) | *S. mutans* vs pond water, Burger model, parallel disks | instantaneous shear modulus 63.9 Pa (*S. mutans*), 1.2 Pa (pond water) | https://scholarworks.montana.edu/items/64e6ffe8-dc82-4137-8bc1-87c3c1514774/full |
| Southampton thesis (eprints 359755) | *S. mutans* | E = 380 Pa (± 350, n = 30), stress relaxation 12 s (± 11) | E rises with sucrose |
| *S. mutans* biofilm AFM (2014, *Biofouling*) | cell wall, Z-ring, EPS | moduli by region | numbers in the full text |

Range: 1e1 to 1e4 Pa depending on method; the 10 Pa of Klempt 2024 lies at the low end. With E = 1 to 10 kPa the model
stresses would be 1e-2 to 1e-1 Pa.

## Growth-induced stress (not oral, but direct evidence)
| source | what | value |
|---|---|---|
| Trejo et al. 2013, PNAS (*B. subtilis* pellicles) | internal compressive stress from growth, wrinkling | about 80 Pa; wrinkle wavelength 1-2 mm |
| Geisel et al. 2022, eLife (*P. aeruginosa*, microfluidic) | growth under lateral confinement -> buckling-delamination | stress inferred from buckling |
| Fei et al. 2019, arXiv:1910.14645 | nutrient depletion -> anisotropic stress, wrinkle pattern | model + agar experiments |

The bars example (thesis Ch. 4: confinement raises the seed stress tenfold) is the analogue of the confinement in these
experiments.

## Antibiotic response of multispecies oral biofilms
| source | what |
|---|---|
| Gränicher et al. 2021, Front. Microbiol. 12:741863 | six-species supragingival (Zurich model), low-dose chlorhexidine, CLSM live/dead after 64 h |
| Hope et al. 2002, Antimicrob. Agents Chemother. | chlorhexidine and oral biofilm vitality over time (viability profiling) |
| Blanc et al. 2014, J. Periodontal Res. | flow-system multispecies model, live/dead CLSM, chlorhexidine rinses |

Candidates for a test of the point model with antibiotic (live fraction over time) once the full texts are read.

## To ask on 1 Dec (Prof. Junker, Dr. Soleimani)
Whether the MHH partners (Heine et al.) have stiffness or thickness data of the same biofilms as the calibration data.
