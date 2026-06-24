"""Anatomically faithful peri-implant TRANSMUCOSAL model (axisymmetric CAX4): titanium implant+abutment,
alveolar bone, peri-implant mucosa (gingiva cuff), the gingival sulcus, and the BIOFILM as a thin layer
on the titanium surface WITHIN the sulcus -- between the Ti and the gingiva, exactly where peri-implant
biofilm sits.  As peri-implantitis advances the pocket deepens and the biofilm DESCENDS apically toward
the bone crest; we parameterise the pocket bottom to show that descent.

argv: pocket_bottom_mm  tag      (z of pocket floor relative to bone crest z=0; +=supracrestal/healthy,
                                   -=infrabony/disease)
Structured (r,z) grid -> CAX4. Units mm, MPa.  Writes tm_<tag>.inp.
"""
import sys
import numpy as np

OUT = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real"
PB = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0      # pocket bottom (mm, rel. bone crest)
TAG = sys.argv[2] if len(sys.argv) > 2 else "p0"

R_I = 2.05            # implant/abutment radius
W_SULC = 0.35         # sulcus width (biofilm layer thickness on Ti)
R_B = 5.0            # outer bone / gingiva radius
Z_BONE, H_AB = -8.0, 4.0    # bone bottom z, abutment top z
H_G = 3.0            # gingiva height above crest (~biologic width)
EPS_GROWTH = 0.19
DR, DZ = 0.12, 0.15
MATS = {"TI": (110000., 0.34), "BONE": (13700., 0.30), "GINGIVA": (3.0, 0.45),
        "BIOFILM": (1.0, 0.45)}


def region(rc, zc):
    if rc < R_I:
        return "TI"                                   # implant (z<0) + abutment (z>0)
    # biofilm fills the peri-implant gap on the Ti surface from the gingival margin (H_G) down to the
    # pocket bottom PB. For PB<0 the pocket is INFRABONY: the bone has resorbed and the defect is
    # colonised by biofilm -> the biofilm DESCENDS below the crest along the implant.
    if R_I <= rc <= R_I + W_SULC and PB <= zc <= H_G:
        return "BIOFILM"
    if zc < 0.0 and rc <= R_B:
        return "BONE"
    if 0.0 <= zc <= H_G and rc <= R_B:
        return "GINGIVA"                              # mucosa cuff outside the sulcus
    return None                                       # air above gingiva


def main():
    Nr = int(round(R_B / DR)); Nz = int(round((H_AB - Z_BONE) / DZ))
    rs = np.linspace(0, R_B, Nr + 1); zs = np.linspace(Z_BONE, H_AB, Nz + 1)

    def nid(i, j):
        return j * (Nr + 1) + i + 1

    elems = {m: [] for m in MATS}
    e = 0
    for j in range(Nz):
        for i in range(Nr):
            rc = 0.5 * (rs[i] + rs[i + 1]); zc = 0.5 * (zs[j] + zs[j + 1])
            reg = region(rc, zc)
            if reg is None:
                continue
            e += 1
            elems[reg].append((e, [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]))

    L = []; ap = L.append
    ap("*HEADING"); ap(" peri-implant transmucosal axisymmetric %s pocket_bottom=%g" % (TAG, PB))
    ap("*NODE")
    for j in range(Nz + 1):
        for i in range(Nr + 1):
            ap(" %d, %.5f, %.5f" % (nid(i, j), rs[i], zs[j]))
    for m, els in elems.items():
        if not els:
            continue
        ap("*ELEMENT, TYPE=CAX4, ELSET=%s" % m)
        for eid, c in els:
            ap(" %d, %d, %d, %d, %d" % (eid, c[0], c[1], c[2], c[3]))
    for m, (E, nu) in MATS.items():
        if not elems[m]:
            continue
        ap("*SOLID SECTION, ELSET=%s, MATERIAL=%s" % (m, m))
        ap("*MATERIAL, NAME=%s" % m); ap("*ELASTIC"); ap(" %.1f, %.3f" % (E, nu))
        if m == "BIOFILM":
            ap("*EXPANSION"); ap(" %.6f" % EPS_GROWTH)

    fixed, abut_top, alln = [], [], []
    for j in range(Nz + 1):
        for i in range(Nr + 1):
            r, z = rs[i], zs[j]; nd = nid(i, j)
            alln.append(nd)
            if z <= Z_BONE + 1e-6 or r >= R_B - 1e-6:
                fixed.append(nd)
            if abs(z - H_AB) < 1e-6 and r < R_I:
                abut_top.append(nd)
    for nm, ids in (("FIXED", fixed), ("ABUT", abut_top), ("ALLN", alln)):
        ap("*NSET, NSET=%s" % nm); ids = sorted(set(ids))
        for k in range(0, len(ids), 16):
            ap(" " + ",".join(str(x) for x in ids[k:k + 16]))

    ap("*INITIAL CONDITIONS, TYPE=TEMPERATURE"); ap(" ALLN, 0.0")
    ap("*BOUNDARY"); ap(" FIXED, 1, 2")
    ap("*STEP"); ap(" 1) biofilm growth"); ap("*STATIC")
    ap("*TEMPERATURE"); ap(" ALLN, 1.0")
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    ap("*STEP"); ap(" 2) occlusal load"); ap("*STATIC")
    ap("*CLOAD")
    nab = max(1, len(set(abut_top)))
    for n in sorted(set(abut_top)):
        ap(" %d, 2, %.4f" % (n, -100.0 / nab))
    ap("*OUTPUT, FIELD"); ap("*NODE OUTPUT"); ap(" U")
    ap("*ELEMENT OUTPUT, POSITION=CENTROID"); ap(" S, COORD"); ap("*END STEP")
    open(f"{OUT}/tm_{TAG}.inp", "w").write("\n".join(L) + "\n")
    tot = sum(len(v) for v in elems.values())
    print("wrote tm_%s.inp (%d CAX4: %s) pocket_bottom=%g"
          % (TAG, tot, {m: len(v) for m, v in elems.items() if v}, PB))


if __name__ == "__main__":
    main()
