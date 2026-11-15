# Chlorhexidine on a six-species oral biofilm (Gränicher et al. 2021) and the antibiotic term of the point model

Note of 13 Nov 2026 for item 4 of "実データ系" in `TODO_2026-11.ja.md`. Sources read: the full text of Gränicher et al.
saved in this folder (`Granicher2021_CHX_multispecies_FrontMicrobiol.txt/.pdf`), and
`Klempt2026_multispecies_interaction_AAM.txt`. Numbers marked (unverified) were not read by me in the source itself.

## 1. Reference (verified against the saved full text)

K. A. Gränicher, L. Karygianni, T. Attin, T. Thurnheer, "Low Concentrations of Chlorhexidine Inhibit the Formation and
Structural Integrity of Enzyme-Treated Multispecies Oral Biofilms", Front. Microbiol. 12:741863 (28 Sep 2021),
doi:10.3389/fmicb.2021.741863, PMCID PMC8506149, CC BY.

## 2. What was done and measured

| item | content | where in the paper |
|---|---|---|
| biofilm | Zurich supragingival model, six species: *S. mutans* OMZ 918, *S. oralis* OMZ 607, *A. oris* OMZ 745, *V. dispar* OMZ 493, *F. nucleatum* OMZ 598, *C. albicans* OMZ 110 (a yeast); pellicle-coated hydroxyapatite disks (9 mm), anaerobic, 37 °C, 64 h | Methods |
| medium | 70 % saliva, 30 % mFUM; 0.3 % glucose (0 to 16 h), then 0.15 % glucose and 0.15 % sucrose; changed at 16 and 40 h | Methods |
| chlorhexidine (CHX) | 0.05, 0.1, 0.2 % (0.2 % = positive control), untreated control; each with or without DNase I (0.002 mg/ml) and proteinase K (0.1 mg/ml) in the medium | Table 1 |
| exposure | six dips of 1 min, at 16, 20, 24, 40, 44, 48 h, each followed by three saline dips; harvest at 64 h (16 h after the last dip) | Methods, Table 1 footnote |
| measured | CFU of the whole biofilm and of each species on selective agar, n = 6; CLSM images of matrix stains and of Live/Dead stain (representative images only, no quantified live fraction, no thickness) | Fig. 1, Table 2, Figs. 2 and 3 |

The exposure covers the growth phase, so the end point is the biofilm that formed under repeated short doses
(the title says "inhibit the formation"), not the kill of a mature biofilm by one dose.

Results in log10 CFU, mean ± SD. Columns 2 and 3 are Table 2 of the paper. Column 1 is not in this paper as numbers
(only box plots, Fig. 1); the values are from the text of the companion paper with the same protocol, Karygianni et al.
2020, J. Clin. Med. 9:983, doi:10.3390/jcm9040983, which I read only through an automated summary of the PMC page
(PMC7231231): (unverified). Column 4 is my subtraction, column 1 minus column 2, so it is a comparison across two papers.

| species | untreated control (unverified) | 0.05 % CHX | enzymes + 0.05 % CHX | reduction by 0.05 % CHX (derived, unverified) |
|---|---|---|---|---|
| total | 7.95 ± 0.16 | 6.44 ± 0.93 | 3.88 ± 1.09 | 1.5 |
| *A. oris* | 6.99 ± 0.36 | 6.01 ± 1.01 | 3.35 ± 0.96 | 1.0 |
| *V. dispar* | not given in text | 5.69 ± 1.35 | 3.46 ± 1.23 | unknown |
| *F. nucleatum* | not given in text | 3.51 ± 1.93 | 2.58 ± 1.11 (p = 0.208) | unknown, the largest of the six |
| *S. mutans* | 6.62 ± 0.23 | 5.05 ± 2.78 | 2.74 ± 0.97 | 1.6 |
| *S. oralis* | 6.94 ± 0.55 (6.84 ± 0.87 elsewhere in the same paper) | 5.44 ± 1.19 | 3.39 ± 1.02 | 1.5 (1.4) |
| *C. albicans* | 3.25 ± 0.59 | 1.98 ± 1.25 | 0.22 ± 0.53 | 1.3 |

Further results and caveats:
- 0.1 % and 0.2 % CHX: no growth detected, all counts below the detection level (Results, text; not shown in Fig. 1).
  So only one dose (0.05 %) gives finite numbers, and the dose response is one point plus a bracket (the step lies
  between 0.05 and 0.1 %).
- The Results text gives a control of 7.95 ± 0.17 log10 CFU in the sentence on *F. nucleatum*. This equals the total
  count of the control in Karygianni 2020 and cannot be *F. nucleatum* alone (the other species already sum to about a
  quarter of that total), so I read it as the total. The *F. nucleatum* control therefore has to be read from Fig. 1.
- Units: "log10 CFU", from one disk harvested into 1 ml saline. Missing values were set to the detection limit.
- The paper writes "0.5 % CHX" in the descriptions of Figs. 2 and 3 and once in the Discussion; from Table 1 this is 0.05 %.
- The SDs are large (up to ± 2.78 for *S. mutans*), with n = 6.
- CLSM (Figs. 2, 3): less dense biofilms after 0.05 % CHX, more so with enzymes; after 0.2 % CHX a few adherent cells
  with damaged membranes remain. Qualitative only.

## 3. How the antibiotic enters the point model (Klempt et al. 2026, Sect. 2)

- State per species: volume fraction $\phi_i$ (living and dead cells), fraction of living cells $\psi_i$, living
  volume fraction $\bar\phi_i=\phi_i\psi_i$ (Eq. 8); dead cells $\phi_i(1-\psi_i)$; empty space $\phi_0$ with
  $\sum_{l=0}^n\phi_l=1$ (Eq. 9). Nutrient $c^*$ and antibiotic $\alpha^*$ are given functions of time.
- Energy (Eq. 10): $-\tfrac12 c^*\,\bar{\boldsymbol\phi}\cdot\mathbf A\cdot\bar{\boldsymbol\phi}
  +\tfrac12\alpha^*\,\boldsymbol\psi\cdot\mathbf B\cdot\boldsymbol\psi$ with the diagonal antibiotic sensitivity matrix
  $\mathbf B=\mathrm{diag}(b_1,\dots,b_n)$ (Eq. 12).
- The antibiotic appears in one place of the evolution equations, the equation for $\psi_i$ (Eq. 17 as printed):
  $0=-c^*\phi_i\big(a_{ii}\bar\phi_i+\sum_j a_{ij}\bar\phi_j\big)+\alpha^*\psi_i b_i
  +\eta_i(\dot\psi_i\phi_i^2+\bar\phi_i\dot\phi_i)+\gamma$. The term $\alpha^*b_i\psi_i$ lowers $\psi_i$; it is linear
  in $\alpha^*$, in $b_i$ and in $\psi_i$, and it is opposed by the nutrient term. The equation for $\phi_i$ (Eq. 16)
  has no antibiotic term.
- Species-specific antibiotic parameter: $b_i$ only. $\eta_i$ is the viscosity (rate of response of species $i$ to
  nutrient and antibiotic alike), not an antibiotic parameter. There is no cross-sensitivity ($\mathbf B$ is diagonal)
  and no threshold or saturation constant.
- What falls: $\psi_i$ directly, hence $\bar\phi_i$. $\phi_i$ changes only indirectly: its growth term in Eq. 16
  carries $\psi_i$, so a killed species stops growing and loses space to others through the constraint. Dead cells stay
  in $\phi_i$; nothing removes them.
- Published cases: $\alpha^*=10$ constant with $b=(1,2)$ (Table 1, cases 4 and 5) and $b=(0.4,0.3,0.2,0.1)$ (Table 2,
  cases 1 to 3); case 4 of Table 2 with $b=(10,2,1,0.01)$ and $\alpha^*=0$, then 100 after one third of the time: species
  1 to 3 die almost at once, species 4 takes over the space, species 2 regains some space (Fig. 11).

## 4. Comparison

What the term can represent as it is:
- Species-specific sensitivity. CFU counts culturable cells, so the survival ratio of species $i$ (treated over control)
  corresponds to the ratio of $\bar\phi_i$ at the end time, if the biofilm volume per disk is taken as equal in both
  groups or the total count is used to scale it. The order at 0.05 % (*A. oris* least reduced, the two streptococci
  about 1.5 log, *F. nucleatum* most) can be set through the order of the $b_i$.
- Dose dependence through $\alpha^*$, and a pulsed schedule through $\alpha^*(t)$, which the paper allows as any function.
- Shift of composition after treatment (the less sensitive species gain share), as in case 4.
- A steep response: my 3D sweep showed the effect switching on between $\alpha^*_0=30$ and 100
  (`PARAMETER_TABLES.md`); the data switch from about 1.5 log to no growth between 0.05 and 0.1 %. This is a
  qualitative parallel only.

What it cannot represent, or what the data do not determine:
- The enzyme effect (2 to 3 log more reduction with DNase I and proteinase K). The model has no matrix and no
  penetration limit; the paper explains this effect by easier diffusion of CHX, without measuring it.
- Loss of biomass and density seen in the CLSM images: in the model dead cells remain in $\phi_i$.
- Survival below about $10^{-2}$ to $10^{-5}$ of the control (0.1 %, 0.2 %): this needs $\psi_i$ resolved very close to
  0, where the barrier on $\psi_i\in[0,1]$ acts, so I would only use "all species dead" as a qualitative check.
- $b_i$ cannot be identified alone. Survival at 64 h depends on $\alpha^*b_i$ together with $c^*$, $\mathbf A$,
  $\eta_i$ and the time scale, and the map from % CHX to $\alpha^*$ is unknown. One finite dose and one end time give
  six numbers for more unknowns than that.
- No time course and no quantified live fraction, so neither $\psi_i(t)$ nor the regrowth between dips is observed.
- *C. albicans* is a yeast and *S. mutans*, *A. oris* are not in the calibrated five-species set (*S. oralis*,
  *A. naeslundii*, *Veillonella*, *F. nucleatum*, *P. gingivalis*); three genera overlap.

Minimal use of these data: fix $\alpha^*_0$ for "0.05 %" by assumption (for example the 100 of case 4), keep
$\mathbf A$, $\eta_i$, $c^*$ of a chosen case, apply $\alpha^*(t)$ as six pulses or as a constant effective level over
16 to 48 h of 64 h (an assumption either way), and choose the ratios $b_i/b_{A.\,oris}$ so that the ratios of
$\bar\phi_i$(treated)/$\bar\phi_i$(control) at the end reproduce the order and roughly the size of the log reductions.
Then use the 0.1 % result (twice $\alpha^*_0$, all species dead) as a check that was not fitted. This is an
illustration of how $b_i$ could be set, not a calibration.

## 5. Other papers (pointers; read only through search results and automated page summaries, so all details unverified)

- Schwendener, Flury, Jenzer, Thurnheer, Karygianni, "PMA-qPCR to quantify viable cells in multispecies oral biofilm
  after disinfectant treatments", Biofilm 9:100281 (2025), doi:10.1016/j.bioflm.2025.100281, PMC12051517. Same Zurich
  model with five bacterial species, 0.2 % CHX for 2 min, once at 64 h or six times; per species CFU, total qPCR and
  PMA-qPCR (viable cells). Total and viable counts per species are the closest measured counterpart of $\phi_i$ and
  $\bar\phi_i$, so this looks better suited than Gränicher. Per-species values are in Figs. 1 and 2, not in tables.
- Chatzigiannidou, Teughels, Van de Wiele, Boon, "Oral biofilms exposure to chlorhexidine results in altered microbial
  composition and metabolic profile", npj Biofilms Microbiomes 6:13 (2020), doi:10.1038/s41522-020-0124-3. Defined
  14-species biofilm, 0.12 % CHX for 5 min every 24 h for three days; intact cell counts over time (about 3 log drop,
  regrowth by 72 h) and strain-level composition (shift from *Veillonella* to *S. gordonii*). A time course with
  regrowth and takeover, the behaviour of case 4.
- Shen et al., Sci. Rep. 6:27537 (2016), doi:10.1038/srep27537: plaque-derived biofilm, 2 % CHX for 1, 3, 10 min, live
  fraction and thickness over eight weeks, with its own ODE model. Not species-resolved.
- Hope and Wilson, Antimicrob. Agents Chemother. 48:1461 (2004), doi:10.1128/AAC.48.5.1461-1468.2004 (listed as "Hope et
  al. 2002" in `EXPERIMENTAL_DATA.md`; the year there should be checked), and Blanc et al., J. Periodontal Res. 49:323
  (2014), doi:10.1111/jre.12110: live/dead profiles, as far as I could see not per species.

## 6. What could be done next

1. Read the per-species control values off Fig. 1 of Gränicher (or ask the authors; the data availability statement
   offers the raw data) and confirm the Karygianni 2020 control numbers in the full text.
2. Read Schwendener et al. 2025 in full and extract viable and total counts per species from Figs. 1 and 2.
3. Run the point model alone with five or six species and pulsed $\alpha^*(t)$, and report which ratios of $b_i$
   reproduce the order of the reductions at one dose; state the assumptions of Sect. 4.
4. Use the 0.1 % result and the composition shift of Chatzigiannidou et al. as qualitative checks.
5. Keep this as a qualitative comparison in the Keio continuation; the model equations stay as published.
