# Geometry files

## `mendeley_xjsx7nfhj8_D46.stl` — lower right first molar (FDI 46), surface mesh

Source: Cristian Diaz, Carlos Andres Ferro Sanchez, Oscar Campo (Universidad
Autonoma de Occidente), *Data of synthetic 3D models of the human jaw,
including teeth, ligaments, and bone structures*, Mendeley Data, V1 (12 Sep 2024),
doi:[10.17632/xjsx7nfhj8.1](https://doi.org/10.17632/xjsx7nfhj8.1).
File `Human Lower Jaw Dataset/Teeth/STL/D46.stl` of the archive
`Human Lower Jaw Dataset.zip` (sha256 of the zip
`902e50c488e8c09bda23bfa6f1d01d5e0ed5528e2738731c2199270ad232e87a`).

Licence: **Creative Commons Attribution 4.0 International (CC BY 4.0)**, as
stated on the dataset page. The file is bundled verbatim (binary STL,
sha256 `652e2b364a3b17334d1f1bfd5769600768cee99322ae52febab257a72a2ae809`);
no change has been made to it. The scripts that use it translate it in memory
only (`ansys_usermat/apdl/anim_real_tooth_3d.py`).

What it is: a closed, manifold triangle surface (5730 points, 11456 faces,
units mm, about 11 x 16 x 19 mm including the two roots; the crown is about
10 x 12 mm). The dataset authors segmented the roots from a CBCT (DICOM) scan
of one adult and replaced the crowns by an intraoral scan (iTero Element 2),
merged in Meshmixer and remeshed in 3-matic; the z axis is the vertical axis
of the jaw. The dataset page gives no statement on consent or anonymisation;
a single tooth surface does not identify a person, and the authors released
the data under CC BY 4.0 themselves.

The other 13 teeth, the periodontal ligaments and the bone of the dataset
are not bundled; download the archive from the DOI if they are needed.
