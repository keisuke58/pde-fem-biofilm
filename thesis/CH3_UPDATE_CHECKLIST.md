# Chapter 3 update checklist (raw-data runs)

Written 16 Oct 2026. Use when the raw-data (measured, two-channel likelihood)
runs of the EML repository (Tmcmc202601) are final. Decision of 16 Oct: DS from
its final stage once it passes; **DH final run = p2** (all criteria pass); the
DH "ult" run (20,000 particles) splits across seeds into different modes and is
reported as a limitation with its numbers (median shifts up to 0.79 posterior sd;
Pg D21/15 = 11.5 / 13.7 / 18.9 against 14.9 measured). Before using DH p2, check
that it uses the same likelihood as ult and that its Pg value agrees across seeds
(`KEIO_PLAN.ja.md` §1, P0).

The same numbers go into the P0 manuscript (Extreme Mechanics Letters, to the
co-authors in the week of 26 Oct), so update both from one source.

## 1. Inputs to collect from the EML repository

- final MAP and posterior samples per condition (CS, CH, DS, DH), commit sha
- criteria (a)-(e) per stage and seed, particles N_p, wall-clock times
- composition RMSE (stage (i) and final), max log-likelihood, ln Z
- identified entries per condition, posterior intervals of A
- pH check (R^2, RMSE, offset in CH)
- figures `eml_*.pdf` regenerated from the final runs

## 2. Places in the thesis to change

`chapters/ch3_heine.tex`
- [ ] l. 246-249: remove the italic note "numbers and figures are from the runs on
      values read off the published figures"
- [ ] l. 252-265: which stages meet (a)-(e); DS/DH status; particle counts
      (18,000 DS, 14,000 DH); composition error ranges; add the DH ult limitation
- [ ] Table `tab:heine_timing` (l. 203-220): N_p ranges and times
- [ ] Table `tab:heine_rmse` (l. 281-299) and its caption (l. 296-297, "provisional")
- [ ] Figure `eml_fit.pdf` (l. 270) and caption l. 277 ("DS and DH from stage (iii)")
- [ ] l. 319-324: number of identified entries (now 4, 6, 8, 9)
- [ ] Figure `eml_ident.pdf` (l. 330), caption l. 337
- [ ] Table `tab:heine_ident` (l. 370-399), caption l. 397-398
- [ ] Figure `eml_posteriors.pdf` (l. 404), caption l. 406-407
- [ ] l. 361-366: signs and intervals of the DS/DH entries
- [ ] Figure `eml_seed_maps.pdf` (l. 421)
- [ ] l. 432-442: ln B values (+0.10, +0.05, -2.67, +1.94) and the Fn-Pg coefficient
- [ ] l. 450-460: pH check (R^2 = 0.63, RMSE 0.12, CH offset 0.13), drop "provisional"
      if DH is final
- [ ] Figure `eml_checks.pdf` (l. 476) and the Pg ratio text (l. 478)
- [ ] l. 47-58: data description (medians vs replicate-level data; live fraction)

Other chapters
- [ ] `main.tex` abstract: counts of identified entries, seeds, speed-up
- [ ] `ch6_conclusion.tex` l. 16-19: DS/DH "provisional", 4/6/8/9 entries
- [ ] `ch6_conclusion.tex` l. 84-86 and l. 125: pH wording (R^2, provisional)
- [ ] `ch2_theory.tex` l. 375: particle counts (5,000 to 18,000)
- [ ] `appendix.tex` (TMCMC settings): particle counts
- [ ] `ch4_ansys.tex` five-species section (source note l. 936): uses the raw CH MAP
      6e18e0b. If the final CH MAP differs, rerun the 8^3 five-species run on
      IKMHIWI03 (`run_calibrated.ps1`) and update the numbers (seed alpha - 1 =
      6.752e-4, vM 4.68e-4 Pa, shares 0.74 -> 0.54 etc.)

## 3. After the update

- [ ] `python3 scripts/audit_thesis_numbers.py` (0 FAIL; Ch. 3 checks only internal
      consistency, the EML repo is not in this checkout)
- [ ] pdflatex + bibtex + pdflatex x2, no errors, page count
- [ ] hand Chapter 3 to Assoc. Prof. Muramatsu (she has Chapters 1 and 2)
