"""extract_klempt_odb.py
=======================
Abaqus Python script: extract element-level SDV + stress from Klempt FEM ODBs.
Run with:  abaqus python extract_klempt_odb.py

Outputs one CSV per (geometry, condition):
  klempt_extract_{geom}_{cond}.csv
  columns: elem_label, phi_gate, alpha, E_gated_MPa, sigma_mises_MPa

phi_gate (SDV5) is used as depth proxy:
  phi_gate = 1.0  ->  inner face (Ti / enamel surface, max growth)
  phi_gate = 0.0  ->  outer face (planktonic, zero growth)

Mode A SDV mapping:
  SDV1=s  SDV2=Je  SDV3=alpha  SDV4=E_gated[MPa]  SDV5=phi_gate
"""

import os
import csv
from odbAccess import openOdb

HERE = os.path.dirname(os.path.abspath(__file__))

CONDS = ["commensal_static", "commensal_hobic",
         "dysbiotic_static", "dysbiotic_hobic"]

# ODB paths: tooth (p23_klempt_A_) and implant (p23imp_klempt_A_)
GEOMS = {
    "tooth":   os.path.join(HERE, "p23_klempt_A_%s.odb"),
    "implant": os.path.join(HERE, "p23imp_klempt_A_%s.odb"),
}

# SDV index (1-based) -> column name
SDV_MAP = {
    "alpha":      3,   # SDV3
    "E_gated":    4,   # SDV4
    "phi_gate":   5,   # SDV5
}


def extract_one(odb_path, out_csv):
    print("Opening: %s" % odb_path)
    odb = openOdb(odb_path, readOnly=True)

    step  = odb.steps.values()[-1]
    frame = step.frames[-1]
    fo    = frame.fieldOutputs

    # build lookup: elemLabel -> index for SDV fields
    fo_s = fo["S"]

    # SDV fields
    fo_sdv = {}
    for name, idx in SDV_MAP.items():
        key = "SDV%d" % idx
        if key not in fo.keys():
            print("  WARNING: %s not found in frame outputs" % key)
            fo_sdv[name] = None
        else:
            fo_sdv[name] = fo[key]

    n_elem = len(fo_s.values)
    print("  Elements: %d" % n_elem)

    rows = []
    for i in range(n_elem):
        sv = fo_s.values[i]
        elem_label = sv.elementLabel
        mises      = sv.mises

        phi_gate = fo_sdv["phi_gate"].values[i].data  if fo_sdv["phi_gate"]  else 0.0
        alpha    = fo_sdv["alpha"].values[i].data     if fo_sdv["alpha"]     else 0.0
        e_gated  = fo_sdv["E_gated"].values[i].data  if fo_sdv["E_gated"]   else 0.0

        rows.append((elem_label, phi_gate, alpha, e_gated, mises))

    odb.close()

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["elem_label", "phi_gate", "alpha", "E_gated_MPa", "sigma_mises_MPa"])
        for r in rows:
            w.writerow(r)

    print("  Wrote %d rows -> %s" % (len(rows), out_csv))


def main():
    for geom, path_template in GEOMS.items():
        for cond in CONDS:
            odb_path = path_template % cond
            if not os.path.exists(odb_path):
                print("SKIP (not found): %s" % odb_path)
                continue
            out_csv = os.path.join(HERE, "klempt_extract_%s_%s.csv" % (geom, cond))
            extract_one(odb_path, out_csv)

    print("\nDone. Run plot_klempt_depth.py for figures.")


if __name__ == "__main__":
    main()
