# References worth adding, and why

Taken from the reference list of Klempt et al. (2024) — the authors' own
selection of the work their model sits on — checked against what this
repository already cites. Grep counts below are from 2026-09-30.

Present already: Rodriguez (4 mentions), Junker (18), Rath (44), Toyofuku (1),
Billings (9), Flemming (2), Albero (1).

---

## 1. The closest existing work, and it is not cited at all

> **Soleimani M, Szafranski SP, Qu T, Mukherjee R, Stiesch M, Wriggers P,
> Junker P (2023).** Numerical and experimental investigation of multi-species
> bacterial co-aggregation. *Scientific Reports* **13**:11839.

**Zero mentions in this repository.** It is *multi-species*, it is *numerical
and experimental*, and the author list is the supervisor, the dental clinic
side and both Klempt co-authors. If any single paper is the precedent this
thesis is a "miniature version" of, it is this one, and a reader from that
group will notice its absence before anything else.

## 2. Growth kinematics — the foundations are missing entirely

The thesis is built on `F = Fe·Fg` with `Fg = (1+α)I`, and cites the Klempt
papers for it. The papers themselves cite where the decomposition comes from,
and none of those are here.

| | mentions |
|---|---|
| **Rodriguez EK, Hoger A, McCulloch AD (1994).** Stress-dependent finite growth in soft elastic tissues. *J Biomech* **27**(4):455–467 | 4 — present, but this is *the* origin of the multiplicative split and should be the citation at the point `Fg` is introduced |
| **Goriely A (2017).** *The Mathematics and Mechanics of Biological Growth.* Springer | **0** |
| **Lubarda VA, Hoger A (2002).** On the mechanics of solids with a growing mass. *Int J Solids Struct* **39**(18):4627–4664 | **0** |
| **Epstein M, Maugin GA (2000).** Thermomechanics of volumetric growth in uniform bodies. *Int J Plast* **16**(7–8):951–978 | **0** |

## 3. The supervisor's own growth FEM line

Both are directly load-bearing for choices this thesis has already made.

- **Soleimani M (2019).** Finite strain visco-elastic growth driven by nutrient
  diffusion: theory, FEM implementation and an application to the biofilm
  growth. *Comput Mech* **64**(5):1289–1301.
  → This repository has a viscoelastic UMAT (`F = Fe·Fv·Fg`, continuation item
  D). This is its precedent, by the supervisor.
- **Soleimani M, Muthyala N, Marino M, Wriggers P (2020).** A novel
  stress-induced anisotropic growth model driven by nutrient diffusion. *J Mech
  Phys Solids* **144**:104097.
  → Stress → growth feedback is exactly what the one-way coupling in the Oliver
  deck's proposal gives up. Citing the paper that models it is the strongest
  way to say what is being deferred and why it is deferrable.

## 4. Evidence that growth-induced stress matters biologically

The thesis computes stress inside a biofilm; the obvious question is what it is
for. This is the paper that answers it, and it is not cited.

- **Chu EK, Kilic O, Cho H, Groisman A, Levchenko A (2018).** Self-induced
  mechanical stress can trigger biofilm formation in uropathogenic
  *Escherichia coli*. *Nat Commun* **9**:4087. — **0 mentions.**

## 5. Oral, experimental, and from the same institute

- **Rath H, Feng D, Neuweiler I, Stumpp NS, Nackenhorst U, Stiesch M (2017).**
  Biofilm formation by the oral pioneer colonizer *Streptococcus gordonii*: an
  experimental and numerical study. *FEMS Microbiol Ecol* **93**(3):fix010.
  → "Rath" greps 44 times here, but check those are this paper and not another
  match; a pioneer-colonizer streptococcus study is directly comparable to
  *S. oralis* leading the succession panel.

## 6. Biofilm material properties — the E_SPEC limitation

`E_SPEC` is disclosed as an assumption carrying about half the headline ratio.
These are the reviews that say why no better number exists.

- **Böl M, Ehret AE, Bolea Albero A, Hellriegel J, Krull R (2013).** Recent
  advances in mechanical characterisation of biofilm and their significance for
  material modelling. *Crit Rev Biotechnol* **33**(2):145–171.
- **Billings N et al. (2015).** *Rep Prog Phys* **78**(3):036601 — already cited.

## 7. The alternatives this work is not

A continuum multi-species model should say what it is choosing against. The
individual-based simulators are the alternative, and none appear here.

- **Mattei MR et al. (2018).** Continuum and discrete approach in modeling
  biofilm development and structure: a review. *J Math Biol* **76**(4):945–1003
  — the review that frames the choice.
- **Lardon LA et al. (2011).** iDynoMiCS. *Environ Microbiol* **13**(9):2416–2434.
- **Li B et al. (2019).** NUFEB: a massively parallel simulator for
  individual-based modelling of microbial communities. *PLoS Comput Biol*
  **15**(12):e1007125.
- **Naylor J et al. (2017).** Simbiotics. *ACS Synth Biol* **6**(7):1194–1210.

## 8. Theory of the Hamilton approach itself

- **Junker P, Balzani D (2021).** An extended Hamilton principle as unifying
  theory for coupled problems and dissipative microstructure evolution.
  *Continuum Mech Thermodyn* **33**(4):1931–1956.
- **Monod J (1949).** The growth of bacterial cultures. *Annu Rev Microbiol*
  **3**(1):371–394 — the Monod relation the growth function uses.

---

## Suggested order

1. **Soleimani et al. 2023** (co-aggregation) — the missing precedent, read first
2. **Soleimani 2019** and **2020** — the supervisor's growth FEM, and what one-way coupling defers
3. **Rodriguez 1994** — move it to where `Fg` is introduced
4. **Chu et al. 2018** — why the stress is worth computing
5. **Goriely 2017** — the standing reference for growth mechanics, for the theory chapter
6. **Mattei 2018** + one individual-based citation — to position the continuum choice

Nothing here is verified beyond the bibliographic details as printed in Klempt
et al. (2024); the papers themselves have not been read in this repository.

---

## Status, 2026-09-30: all of the above are now in the bibliography

Fourteen entries added to `biofilm_3tooth_refs.bib` (23 → 37), the file
`thesis_ch5/_build_check.tex` and `biofilm_3tooth_report.tex` both use. Brace
balance, duplicate keys and required fields checked.

**Provenance is recorded per entry, because it matters for a submission.**

| verified first-hand (PDF read, DOI included) | key |
|---|---|
| Soleimani et al. 2023, co-aggregation | `Soleimani2023CoAggregation` |
| Soleimani 2019, viscoelastic growth | `Soleimani2019ViscoelasticGrowth` |
| Soleimani, Haverich & Wriggers 2021 | `Soleimani2021Atherosclerosis` |
| Chu et al. 2018 | `Chu2018SelfInducedStress` |

The other ten are **transcribed from Klempt et al. (2024)'s reference list and
have not been checked against the articles**. They deliberately carry no `doi`
field rather than a guessed one, and each entry's `note` says so:
`Soleimani2020AnisotropicGrowth`, `Goriely2017BiologicalGrowth`,
`LubardaHoger2002GrowingMass`, `EpsteinMaugin2000VolumetricGrowth`,
`Rath2017StreptococcusGordonii`, `Bol2013BiofilmMechanicalCharacterisation`,
`Mattei2018ContinuumDiscreteReview`, `Lardon2011iDynoMiCS`, `Li2019NUFEB`,
`Naylor2017Simbiotics`. **Verify these before submission** — a wrong volume or
page range in a thesis bibliography is the kind of error a viva finds.

### The Rath suspicion was right

§5 above wondered whether the 44 "Rath" matches in this repository were the
2017 *S. gordonii* paper. They were not: **`biofilm_3tooth_refs.bib` had no
Rath entry at all** before today. Whatever those matches were, none of them
was a citation. Now added as `Rath2017StreptococcusGordonii`.

### Where each should be cited

Adding the entries is not citing them. The points that need a `\cite`:

- `Soleimani2021Atherosclerosis` — where `Fg = (1+α)I` is introduced (its
  Eq. 15 is this convention verbatim, `α(0) = 0`), and where the growth cap is
  described (its Eq. 17). This is what closes `CITATION_AUDIT.md` F1c.
- `Rodriguez1994StressDependentGrowth` (already present),
  `LubardaHoger2002GrowingMass`, `Goriely2017BiologicalGrowth` — at the
  multiplicative split itself, which currently cites only the Klempt papers.
- `Soleimani2023CoAggregation` — in the introduction, as the precedent this
  work is a smaller version of. A reader from that group will look for it
  first.
- `Soleimani2019ViscoelasticGrowth` — at the viscoelastic UMAT, and again at
  the integrator discussion (`VISCOUS_UPDATE_SCHEME.md`).
- `Soleimani2020AnisotropicGrowth` — where the one-way coupling is justified,
  to say what is being deferred.
- `Chu2018SelfInducedStress` — where the stress magnitudes are interpreted,
  with the `E_SPEC` caveat from `CHU2018_NOTES.md` attached.
- `Bol2013BiofilmMechanicalCharacterisation` — at the `E_SPEC` limitation.
- `Mattei2018ContinuumDiscreteReview`, `Lardon2011iDynoMiCS`, `Li2019NUFEB`,
  `Naylor2017Simbiotics` — where the continuum approach is chosen, to say what
  it is chosen against.
- `EpsteinMaugin2000VolumetricGrowth` — at the thermodynamics of an open
  growing system.
- `Rath2017StreptococcusGordonii` — beside the *S. oralis* pioneer-colonizer
  discussion.
