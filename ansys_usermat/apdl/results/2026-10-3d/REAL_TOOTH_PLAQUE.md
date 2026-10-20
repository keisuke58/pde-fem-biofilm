# Towards a realistic plaque shape on the real molar: three synthetic steps (20 Oct 2026)

Script `ansys_usermat/apdl/real_tooth_plaque.py`; numbers in `real_tooth_plaque.json` (this directory); videos
`assets/anim_real_tooth_plaque_1_3d.mp4` and `assets/anim_real_tooth_plaque_3_3d.mp4` (last frames `.png`); figure
`assets/fig_real_tooth_plaque.pdf` / `.png`; caches `figdata/real_tooth_plaque_*.npz` (git-ignored, rebuilt by the
script when missing). This is a demonstration on the molar FDI 46 of `REAL_TOOTH.md`, not a thesis result, and it uses
no experimental data of the project beyond aggregate statistics of the CLSM stacks (below).

The idealised film of `REAL_TOOTH.md` (uniform 0.5 mm, one flank, 120 degree sector, cervical third, phi = 1 inside) is
replaced in three steps, each built synthetically: (1) where plaque accumulates on a molar, (2) how its thickness varies
along the surface, (3) what it looks like inside. Each step keeps the growth model and the mechanics of the
demonstration; only the initial phi changes.

## Sources

Published values (the full citations; what each one is used for is said at the parameter):

- Al-Ahmad A., Wunder A., Auschill T. M., Follo M., Braun G., Hellwig E., Arweiler N. B. (2007). The in vivo dynamics of
  Streptococcus spp., Actinomyces naeslundii, Fusobacterium nucleatum and Veillonella spp. in dental plaque biofilm as
  analysed by five-colour multiplex fluorescence in situ hybridization. *Journal of Medical Microbiology* 56(5), 681-687.
  In situ biofilm on bovine enamel slabs worn buccally: thickness 14.9 +- 5.0 um after 1 day, 49.3 +- 11.6 um after
  7 days (abstract).
- Wood S. R., Kirkham J., Marsh P. D., Shore R. C., Nattress B., Robinson C. (2000). Architecture of intact natural
  human plaque biofilms studied by confocal laser scanning microscopy. *Journal of Dental Research* 79(1), 21-27.
  Four-day natural plaque on enamel in situ, eight volunteers: thickness 35-215 um at the centre of the devices and
  75-220 um at the sheltered enamel/nylon junction; heterogeneous structure with channels and fluid-filled voids
  (abstract).
- Rabe P., Twetman S., Kinnby B., Svensater G., Davies J. R. (2015). Effect of fluoride and chlorhexidine digluconate
  mouthrinses on plaque biofilms. *The Open Dentistry Journal* 9, 106-111. Seven-day in vivo biofilms on enamel discs:
  28-96 um thick, open architecture with mushroom-like structures in one subject, biomass largest in the middle layers.
- Furuichi Y., Lindhe J., Ramberg P., Volpe A. R. (1992). Patterns of de novo plaque formation in the human dentition.
  *Journal of Clinical Periodontology* 19(6), 423-433. Plaque accumulated most at the approximal surfaces and least at
  the palatal surfaces; the mandible harboured more plaque than the maxilla (abstract).
- Lie T., Gusberti F. (1979). Replica study of plaque formation on human tooth surfaces. *Acta Odontologica
  Scandinavica* 37(2), 65-72. On buccal surfaces bacterial accumulation occurred near the cemento-enamel junction in
  6-hour specimens and expanded in a coronal direction (abstract).
- Loe H., Theilade E., Jensen S. B. (1965). Experimental gingivitis in man. *Journal of Periodontology* 36(3), 177-187.
  The classical observation that plaque forms along the gingival margin when oral hygiene stops.
- Quirynen M., van Steenberghe D. (1989). Is early plaque growth rate constant with time? *Journal of Clinical
  Periodontology* 16(5), 278-283. Four days without hygiene, buccal surfaces of anterior teeth and premolars: 27.7 % of
  the surface area covered on average; the growth pattern follows the irregularities of the tooth surface (abstract).
- Ekstrand K. R., Bjorndal L. (1997). Structural analyses of plaque and caries in relation to the morphology of the
  groove-fossa system on erupting mandibular third molars. *Caries Research* 31(5), 336-348. Microorganisms occupied
  nearly all grooves; viable organisms mainly at the entrance (abstract). Ekstrand K. R., Carlsen O., Thylstrup A.
  (1991). Morphometric analysis of occlusal groove-fossa-system in mandibular third molar. *Scandinavian Journal of
  Dental Research* 99(3), 196-204: interlobal groove depth 0.13-0.55 of the maximum crown height; a groove with a
  structure angle of 25 degrees or less is a fissure (abstract).
- Lawrence J. R., Korber D. R., Hoyle B. D., Costerton J. W., Caldwell D. E. (1991). Optical sectioning of microbial
  biofilms. *Journal of Bacteriology* 173(20), 6558-6567. Biofilms of Pseudomonas and Vibrio by confocal microscopy:
  73-98 % extracellular material and space; Pseudomonas biofilms densest at the attachment surface and more diffuse
  outward (abstract).
- Heydorn A., Nielsen A. T., Hentzer M., Sternberg C., Givskov M., Ersboll B. K., Molin S. (2000). Quantification of
  biofilm structures by the novel computer program COMSTAT. *Microbiology* 146(10), 2395-2407. The descriptors used
  below (mean thickness, roughness coefficient, substratum coverage, biovolume).
- Banas J. A., Fountain T. L., Mazurkiewicz J. E., Sun K., Vickerman M. M. (2007). Streptococcus mutans glucan-binding
  protein-A affects Streptococcus gordonii biofilm architecture. *FEMS Microbiology Letters* 267(1), 80-88. COMSTAT
  values of monospecies S. gordonii biofilms (Table 2): substratum coverage 21-87 %, roughness coefficient 0.40-1.43,
  the rougher biofilms being the less covering ones; used only to place the CLSM statistics below in a published range.

CLSM/FISH stacks of the in-vitro biofilms (CH, DH), N. Heine et al., processed in `nife_intern`
(`results/fish_3d/ic/phi_xyz_{CH,DH}_d{1,6,10,15,21}.npy`, five-species occupancy per voxel, xy 1.44 um, z 2 um,
field of view 184 um): these are in-vitro 21-day films on titanium, about 40 um thick, so only relative quantities
are transferred to the 0.25-0.5 mm plaque on the molar, never the absolute thickness. Aggregates I derived (occupied =
total occupancy above half the stack's 95th percentile; the alternative threshold of a quarter in brackets; the eight
stacks of day 6 and later; the stacks have no visible substratum plane, the signal attenuates with depth, so the
thickness is the biovolume per column, COMSTAT's mean thickness, and the porosity is taken inside each column's
occupied band):

| statistic | day >= 6, mean +- sd over 8 stacks | used as |
|---|---|---|
| coefficient of variation of the local thickness (3 x 3 median, 4.3 um) | 0.26 +- 0.05 (0.16 +- 0.05) | COV_T = 0.25 (step 2); Al-Ahmad's 7-day sd/mean is 0.24 |
| lateral correlation length of the thickness (1/e of the autocorrelation) | 14 +- 7 um (18 +- 9 um), 0.7-1.1 x the mean column thickness, limited by the 184 um field of view | ELL_REL = 1.0 x mean thickness (step 2) |
| porosity inside the occupied band | 0.46 +- 0.22 (0.36 +- 0.17) | POROSITY = 0.40 (step 3) |
| lateral coverage (columns with biomass) | 0.97 +- 0.06 (1.00) | P_BASE = 0.03, the void fraction of the first cell layer (step 3) |
| correlation length of the column signal | 18 +- 8 um, about 0.45 x the occupied band thickness | ELL_C: 0.45 x 0.32 mm = 0.14 mm is below two cells, so ELL_C = 2 h_f = 0.17 mm (step 3) |

No per-stack value and no array of the stacks is in this repository; the JSON holds the aggregates above.

(Results follow once the runs are complete.)
