# Neighbored Element Method (NEM): the papers, and how they bear on this work

Felix Klempt sent the papers on the Neighbored Element Method (NEM) that he
knows of (mail of 6 Oct 2026). NEM is the discretisation behind the partner
element and behind his own code (`ansys_usermat/FELIX_VS_PARTNER_DIFF.md`).
Four of the seven papers are open access and are bundled here as PDF with a
searchable text extraction; licences are listed in `THIRD_PARTY.md`. The
other three are not open access: only their citation and a summary in my own
words are given here, and the PDFs stay outside the repository
(`F:\felix_private\` on IKMHIWI03).

## The idea in one paragraph

The group's material models (topology optimisation, gradient-enhanced damage,
biofilm growth) each have a field variable with a gradient regularisation:
density, damage or phi, governed by a parabolic equation with a Laplacian.
Instead of making that variable a nodal unknown of the finite element system,
NEM keeps it at the integration points (or element centres) as an internal
variable. Its Laplacian and gradient are built from the values at the
neighbouring elements, by a Taylor expansion fitted in the weighted
least-squares sense. That yields an operator matrix that depends on the mesh
only, so it is set up once. The field is then updated explicitly, staggered
from the displacement solve, so the finite element system keeps only the
displacements and costs about as much as a purely elastic one. The price is a
time step limit of the explicit update and an operator that is a
finite-difference-like approximation rather than a Galerkin one.

## The papers (oldest first)

| paper | in repo | what it adds |
|---|---|---|
| Jantos, Hackl, Junker (2019), *An accurate and fast regularization approach to thermodynamic topology optimization*, Int J Numer Methods Eng 117:991-1017, doi:10.1002/nme.5988 | no (not open access) | The starting point. Topology optimisation from Hamilton's principle gives a parabolic evolution equation for the density; it is regularised by its gradient and integrated explicitly. The Laplacian is discretised by a Taylor series over neighbouring elements, giving a constant operator matrix, which also works on unstructured meshes. The regularisation parameter controls the minimum member size a priori. |
| Junker, Schwarz, Jantos, Hackl (2019), *A fast and robust numerical treatment of a gradient-enhanced model for brittle damage*, Int J Multiscale Comput Eng, doi:10.1615/IntJMultCompEng.2018027813 | no (no open licence) | Carries the same treatment over to gradient-enhanced brittle damage, again from Hamilton's principle for non-conservative continua. Here the combination of finite elements and the neighbour-based update gets the name Neighbored Element Method. Several boundary value problems give mesh-independent results at a cost close to an elastic computation. |
| Vogel, Junker (2020), *Adaptive and highly accurate numerical treatment for a gradient-enhanced brittle damage model*, Int J Numer Methods Eng 121:3108-3131, doi:10.1002/nme.6349 | **yes** (CC BY) | The same damage model solved with space- and time-adaptive finite elements, damage resolved down to crack-like bands. A reference solution and a convergence study rather than NEM itself; useful as the "Galerkin" counterpart of NEM results. |
| Blaszczyk, Jantos, Junker (2022), *Application of Taylor series combined with the weighted least square method to thermodynamic topology optimization*, Comput Methods Appl Mech Eng 393:114698, doi:10.1016/j.cma.2022.114698 | no (Elsevier, all rights reserved) | The general form of the operator: Taylor series plus weighted least squares for values given at arbitrary point clouds, including the Neumann boundary condition. The earlier approaches failed on some mesh types; this one is shown to be mesh-independent across element types. This is the "MLS-based finite-difference operator" of the partner element. |
| Junker, Riesselmann, Balzani (2022), *Efficient and robust numerical treatment of a gradient-enhanced damage model at large deformations*, Int J Numer Methods Eng 123:774-793, doi:10.1002/nme.6876 | **yes** (CC BY-NC) | NEM at finite strains, with element erosion for fully damaged material; snapback and springback examples. Shows the staggered explicit update stays robust at large deformation. |
| Rudolf, Klempt, Kök, Soleimani, Jantos, Junker (2025), *Computational efficiency and accuracy of the Neighbored Element Method*, Finite Elem Anal Des 249:104353, doi:10.1016/j.finel.2025.104353 | **yes** (CC BY, added earlier) | The method paper closest to this work: accuracy and cost of NEM against standard finite elements, by authors of the biofilm element. Reading notes: `KEIO_PLAN.ja.md` §4d. |
| von Zabiensky, Jantos, Junker (2026), *A fast and robust third medium contact approach using the neighbored element method*, Finite Elem Anal Des 255:104489, doi:10.1016/j.finel.2025.104489 | **yes** (CC BY-NC-ND) | Third-medium contact: the regularising gradients become extra unknowns at the quadrature points, discretised NEM-style and solved staggered, so linear elements suffice. Shows NEM beyond damage and density, for gradient regularisation in general. |

## What this means for this work

- **Why ANSYS and Abaqus are not expected to agree to the last digit**
  (`abaqus_composition/README.md`, verification logic): the partner element
  keeps phi at the Gauss points with a NEM operator and an explicit update;
  Abaqus solves phi as a nodal Galerkin field, implicitly. Blaszczyk 2022 and
  Rudolf 2025 give the accuracy of the NEM operator, so the differences seen
  under refinement can be judged against them rather than only against each
  other.
- **Time step**: the explicit NEM update has a stability limit tied to beta
  and the element size (Jantos 2019); the partner's front term adds an
  advective limit on top (`KLEMPT2024_REPRODUCTION.md` sec. 8, cell Péclet
  number).
- **Keio**: Vogel 2020 is the pattern for a Galerkin reference against NEM;
  von Zabiensky 2026 shows how a gradient-type field can be added in NEM style
  if the stress-to-growth coupling (`KEIO_PLAN.ja.md`) needs one.
