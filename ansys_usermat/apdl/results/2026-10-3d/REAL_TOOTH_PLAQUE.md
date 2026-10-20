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

## Orientation of the tooth in the box

z is vertical (occlusal plane horizontal, crown top at z = 5.5 mm, CEJ near z = -1.6 mm, as in `REAL_TOOTH.md`). Which
horizontal axis is mesio-distal I decided from two independent features of the STL: the principal axis of the crown
outline between z = 1 and 3 mm is (x, y) = (-0.18, 0.98), 11 degrees off the y axis, with the crown 11.9 mm long along
it and 9.8 mm across (a lower first molar is about 11 mm mesio-distally and 10.5 mm bucco-lingually); and the two roots
of the lower molar, which lie mesially and distally, separate below z = -5 mm along (-0.35, 0.94), 20 degrees off y.
Both say the mesio-distal axis is close to y; the script takes the crown axis. The root at +y inclines further towards
+y with depth, the distal inclination of lower molar roots, so +y is distal and -y mesial. Buccal and lingual (the
+x and -x flanks; the idealised film sat on -x) I do not tell apart: the model gives both the same thickness. The
groove-fossa system I take from the geometry itself (below), not from the orientation.

## How each step is built

All three steps share the field grid (192^3 cells, h_f = 0.083 mm, the grid of `REAL_TOOTH_SUBMODEL.md`), the growth
model and its constants (RK4, masked Laplacian with zero flux at the walls and at the voxel tooth, k_alpha = 1e-3 /T*,
penalty 5, T* = 1.1) and one diffusion coefficient: beta = 0.02 (t_mean / 0.25 mm)^2 with t_mean = 0.321 mm, the mean
thickness of the step-1 map (volume divided by the covered area), i.e. **beta = 0.033 mm^2/T*** (the rule of
`REAL_TOOTH.md`, an assumption, as there). The field cannot resolve a film thinner than three cells, so every map value
between 0.125 mm (T_CUT) and 0.25 mm (T_MIN) is raised to 0.25 mm and below 0.125 mm there is no film; the fraction of
the band at the floor is reported.

**Step 1, distribution.** The thickness map t(theta, z) is built from the azimuth theta of each cell about the crown
axis and its height z, and the film is the set of open cells with signed distance 0 <= sdf < t:

- a cervical band all round the crown from the gingival margin upwards. The margin is placed 0.5 mm coronal to the CEJ
  (Z_GM = -1.1 mm; in health the margin lies coronal to the CEJ; the exact height is my choice). Plaque starts at the
  gingival margin and the CEJ and spreads coronally (Loe et al. 1965; Lie and Gusberti 1979), so the thickness is largest
  at the margin and decays coronally: t = T_BL exp(-((z - Z_GM) / 1.5 mm)^2) on the buccal and lingual surfaces, which
  with T_CUT ends the band 1.25 mm above the margin (the cervical third of the 7 mm crown). Nothing is placed below the
  margin: the gingiva is not modelled and the lower edge of the film is free.
- thicker on the approximal surfaces: w(theta) = exp(-(dtheta / 25 deg)^2) with dtheta the angular distance to the
  mesial or distal direction (the sector half width is my choice), t = w T_IP f_ip + (1 - w) T_BL f_bl. Approximal
  surfaces harbour the most plaque (Furuichi et al. 1992), and there the film reaches up to the contact area, Z_HOC =
  2.0 mm (the crown is widest at z = 1..3 mm in this STL), with f_ip = 1 between the margin and Z_HOC and a decay of
  0.75 mm above it.
- thickness values: the in situ measurements on enamel are 15-50 um over 1-7 days on buccal slabs (Al-Ahmad et al.
  2007), 35-215 um at the centre and 75-220 um at the sheltered edge after 4 days on natural teeth (Wood et al. 2000),
  28-96 um after 7 days (Rabe et al. 2015). All of these are below the three cells of the grid, so **T_BL = 0.25 mm is
  the floor, not a literature value**; the buccal/lingual band is 1.1 to 7 times thicker than measured plaque on smooth
  enamel. **T_IP = 0.5 mm is my choice**: twice T_BL for the approximal surfaces, which hold the most plaque, and the
  thickness of the idealised film, so that the two remain comparable; no in vivo measurement of approximal plaque
  thickness was found (the only values are from devices on buccal surfaces).
- deposits in the occlusal groove-fossa system: the cells above z = 2 mm that a morphological closing of the voxel
  tooth with a ball of radius 0.6 mm fills, i.e. every occlusal concavity narrower than about 1.2 mm is filled with
  plaque (microorganisms occupy nearly all grooves, Ekstrand and Bjorndal 1997; the groove depth is 0.13-0.55 of the
  crown height, Ekstrand et al. 1991; the radius is my choice). This STL has shallow grooves (a few flat facets, see
  `REAL_TOOTH.md`), so the fill is 1.1 mm^3, mostly the mesial fossa.

**Step 2, thickness variation.** t is multiplied by a lognormal field exp(s g - s^2/2), g a Gaussian random field with
covariance exp(-r^2 / ell^2), ell = 1.0 x t_mean = 0.32 mm, s = sqrt(ln(1 + COV_T^2)) with COV_T = 0.25 (from the CLSM
statistics above; Al-Ahmad's 7-day sd/mean is 0.24), seed 0. g is generated on the 3D grid and evaluated at the foot
point x - sdf grad(sdf) of each cell, so it is a field on the tooth surface and the local thickness is constant across
the film. The floor applies after the multiplication (the fraction clipped is reported).

**Step 3, internal structure.** Inside the step-2 envelope a second Gaussian field (covariance length ELL_C = 2 h_f =
0.17 mm, seed 1) sets phi = 0 where it falls below the quantile of a porosity that rises linearly from 0.03 at the
tooth (the stacks' lateral coverage of 0.97-1.0) to 0.77 at the outer surface, mean 0.40 (the stacks' 0.36-0.46), so
that the base layer is almost closed, channels and voids open towards the outside and the outer surface breaks up into
colonies (the "mushroom" tops; Rabe et al. 2015, Wood et al. 2000, Lawrence et al. 1991 for the denser base). The first
cell layer keeps the base porosity. The stiffness E = 10 Pa (phi^2 + 1e-3) then varies inside the film (zero in the
channels until the diffusion fills them).

**Stress.** The film submodel of `real_tooth_submodel.py` on the final state (rigid tooth, crop box around the film,
voxel wall, pyamg + CG to 1e-10): at h = h_f = 0.083 mm with the diffuse front kept (void phi <= 0.05 dropped; level
r1), at h_f with the front dropped (phi <= 0.3; r1t) and at h = h_f/2 = 0.042 mm with the front dropped (r2t): six
elements across the 0.25 mm band and twelve across the approximal film. The front has to go on the fine mesh because
with it the meshes hold 3.1-3.7 M dof (over 11 GB), without it 1.0-1.7 M; r1 -> r1t gives the effect of the front, r1t
-> r2t the mesh refinement at a fixed element set. Regions (field cells at the final state, every stress element
weighted by its fraction of cells in the set, as in the earlier studies): film = phi > 0.3 and open; cervical band = film
with z in [Z_GM, Z_GM + 1 mm] all round; approximal = film with w > 0.5 (within about 21 degrees of the mesio-distal
axis) outside the fissures; buccal/lingual = the rest outside the fissures; fissures = film in the closing fill; the
0.1 mm next to the tooth. Means and percentiles only, as `REAL_TOOTH_SUBMODEL.md` found the pointwise maximum does not
converge (staircase corners).

(Results follow once the runs are complete.)
