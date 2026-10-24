# Note for Oliver (draft, 24 Oct 2026; not sent)

What I found in the element while checking the coupled runs against two other implementations (a finite-volume code
and a Galerkin user element). Short version first, details below. All of it concerns the numerics of the field, not
the model; the fixes I used are switches in USSFin (patches in `oliver_patches/`), nothing is changed by default.

## Short version

1. **Surface points.** phi at the surface points was never updated, so every face acted as phi = 0. With the
   inverse-square mean of the neighbouring Gauss points the walls still lose biomass; taking the nearest Gauss
   points (zero normal gradient) works.
2. **Laplacian.** The MLS Laplacian has zero row sums (constants are exact) but its column sums are not zero, so the
   explicit update does not conserve phi (leak through the walls). A symmetrised copy with the diagonal from the row
   sums (read from a file) conserves it. For that the neighbours must be chosen by distance, not by element order:
   with the element-order choice the 30-point stencils are asymmetric inside on 16^3 and 24^3.
3. **Edges and corners.** With a manufactured solution (known phi and c, sources added) the element does not converge
   from 8^3 to 16^3: the error sits at the edges and corners of the box, where the Neumann condition of the nutrient
   follows the point's own normal (diagonal at a corner) instead of each face. The interior is fine.
4. **Front term.** With beta = 0.02 mm^2/T* and r = 10 mm/T* the cell Peclet number is about 30 and the central MLS
   gradient oscillates. A constant added diffusion r h/4 keeps it stable but smears the seed edge, and the seed
   stress then comes out at a quarter of the converged value on 16^3. One-sided (upwind) MLS gradients of phi, built
   from the same stencils, give a stress within 25 % of the converged value on 16^3 without added diffusion.

## Numbers (separate example with antibiotic, seed means at T* = 1)

| | alpha - 1 | von Mises [Pa] |
|---|---|---|
| element, r h/4, 16^3 | 5.7e-4 | 1.3e-4 |
| element, upwind gradients, 16^3 | 8.0e-4 | 4.6e-4 |
| converged (two second-order codes, 32^3 to 128^3) | 9.0 to 9.3e-4 | 6.0 to 6.3e-4 |

The growth variable and the composition agree between the codes within their mesh-error bands; only the stress is
sensitive to the stabilisation of the front term.

## What I would suggest

- A fit of the surface value that is exact for linear fields (as for the temperature) instead of the mean of the
  neighbours, and a face-wise zero-flux condition at edges and corners.
- Neighbours by distance (one line in GetNhgb).
- An upwind option for the front term, or a limiter on the added diffusion.

Happy to show the runs and the patches whenever it suits you.
