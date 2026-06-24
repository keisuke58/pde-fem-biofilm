"""compare_tooth_klempt.py
==========================
Extract SDV and stress from tooth Klempt UMAT ODB files.

Mode A SDVs: SDV1=s=1+alpha, SDV2=Je, SDV3=alpha, SDV4=E_gated [MPa], SDV5=phi_gate
Mode C SDVs: SDV1=alpha_total, SDV2=Je, SDV3=E_gated, SDV4=s, SDV5=alpha_So, SDV6=alpha_Vd, SDV7=phi_gate

Run via Abaqus Python:
  abaqus python compare_tooth_klempt.py
"""

from odbAccess import openOdb
import numpy as np
import json, os

CONDITIONS = ["commensal_static", "commensal_hobic", "dysbiotic_static", "dysbiotic_hobic"]
MODES      = ["A", "C"]

HERE = os.path.dirname(os.path.abspath(__file__))

def sdv_idx(mode):
    """Return dict mapping name -> 0-based SDV index."""
    if mode == "A":
        return {"s": 0, "Je": 1, "alpha": 2, "E_gated": 3, "phi_gate": 4}
    else:  # C
        return {"alpha_total": 0, "Je": 1, "E_gated": 2, "s": 3,
                "alpha_So": 4, "alpha_Vd": 5, "phi_gate": 6}

def extract_odb(odb_path, mode):
    """Extract last frame SDV + Mises stress statistics."""
    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps.values()[-1]
    frame = step.frames[-1]
    fo = frame.fieldOutputs

    result = {"mode": mode, "odb": os.path.basename(odb_path)}

    # --- SDVs: stored as SDV1, SDV2, ... individual field outputs ---
    sdv_map = sdv_idx(mode)
    for name, idx in sdv_map.items():
        key = "SDV%d" % (idx + 1)
        if key in fo:
            col = [v.data for v in fo[key].values]
            result[name + "_max"]  = float(max(col))
            result[name + "_mean"] = float(sum(col) / len(col))
            result[name + "_min"]  = float(min(col))

    # --- Mises stress and S11 ---
    if "S" in fo:
        mises_list = [v.mises for v in fo["S"].values]
        s11_list   = [v.data[0] for v in fo["S"].values]
        result["mises_max"]  = float(max(mises_list))
        result["mises_mean"] = float(sum(mises_list) / len(mises_list))
        result["s11_max"]    = float(max(s11_list))
        result["s11_min"]    = float(min(s11_list))

    odb.close()
    return result

summary = {}
for mode in MODES:
    summary[mode] = {}
    print(f"\n{'='*60}")
    print(f"  MODE {mode}")
    print(f"{'='*60}")
    print(f"  {'Condition':<22} {'alpha_max':>10} {'E_gated_max':>12} {'phi_gate_max':>13} {'Mises_max [MPa]':>16}")
    print(f"  {'-'*75}")

    for cond in CONDITIONS:
        odb_path = os.path.join(HERE, f"p23_klempt_{mode}_{cond}.odb")
        if not os.path.exists(odb_path):
            print(f"  {cond:<22}  MISSING")
            continue
        r = extract_odb(odb_path, mode)
        summary[mode][cond] = r

        alpha_key = "alpha_max" if mode == "A" else "alpha_total_max"
        alpha_val = r.get(alpha_key, r.get("alpha_max", float("nan")))
        eg_max    = r.get("E_gated_max", float("nan"))
        pg_max    = r.get("phi_gate_max", float("nan"))
        mises_max = r.get("mises_max", float("nan"))

        print(f"  {cond:<22} {alpha_val:>10.4f} {eg_max:>12.6f} {pg_max:>13.4f} {mises_max:>16.4e}")

# Key ratios
print(f"\n{'='*60}")
print("  Key clinical ratio: α_max (CH / DH)")
for mode in MODES:
    try:
        a_key = "alpha_max" if mode == "A" else "alpha_total_max"
        ch = summary[mode]["commensal_hobic"].get(a_key, summary[mode]["commensal_hobic"].get("alpha_max", None))
        dh = summary[mode]["dysbiotic_hobic"].get(a_key, summary[mode]["dysbiotic_hobic"].get("alpha_max", None))
        if ch and dh:
            print(f"  Mode {mode}: CH={ch:.4f}  DH={dh:.4f}  ratio={ch/dh:.2f}x")
    except Exception:
        pass

print(f"\n  E_gated_max (phi^2-gate verification):")
for mode in MODES:
    for cond in CONDITIONS:
        eg = summary[mode].get(cond, {}).get("E_gated_max", None)
        pg = summary[mode].get(cond, {}).get("phi_gate_max", None)
        if eg is not None:
            print(f"  Mode {mode} {cond:<22}: E_gated={eg:.3e} MPa  phi_gate={pg:.4f}")

# Save JSON
out = os.path.join(HERE, "tooth_klempt_comparison.json")
with open(out, "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n  Saved: {out}")
