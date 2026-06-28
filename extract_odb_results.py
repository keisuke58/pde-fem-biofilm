"""
extract_odb_results.py  — run with: abaqus python extract_odb_results.py
Extract S (all components), Mises, U3 from 2-ch UMAT ODBs.
"""
import sys, os, json

try:
    from odbAccess import openOdb
except ImportError:
    print("Must run with: abaqus python", file=sys.stderr)
    sys.exit(1)

RUN_DIR = os.path.expanduser("~/IKM_Hiwi/FEM/_abaqus_runs")
CASES   = [("visco_2ch_ch", "biofilm_visco_2ch_ch"),
           ("visco_2ch_dh", "biofilm_visco_2ch_dh")]


def to_py(obj):
    if isinstance(obj, dict):
        return {k: to_py(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_py(x) for x in obj]
    try:
        return float(obj)
    except (TypeError, ValueError):
        return obj


def extract(case_dir, job_name):
    odb_path = os.path.join(case_dir, job_name + ".odb")
    odb = openOdb(odb_path, readOnly=True)
    step = odb.steps["GROWTH_STATIC"]

    inst_key = list(odb.rootAssembly.instances.keys())[0]
    inst = odb.rootAssembly.instances[inst_key]

    # element centroid z-coordinates (in mm)
    node_z = {nd.label: nd.coordinates[2] for nd in inst.nodes}
    el_z = {}
    for el in inst.elements:
        el_z[el.label] = sum(node_z[n] for n in el.connectivity) / len(el.connectivity)
    el_sorted = sorted(el_z.keys(), key=lambda k: el_z[k])

    times, S11_hist, S22_hist, S33_hist, SMIS_hist, U3_hist = [], [], [], [], [], []

    for frame in step.frames:
        times.append(frame.frameValue)
        fo_keys = list(frame.fieldOutputs.keys())

        def get_fo(name):
            return frame.fieldOutputs[name] if name in fo_keys else None

        sf = get_fo("S")
        uf = get_fo("U")

        s11 = {v.elementLabel: v.data[0] for v in (sf.values if sf else [])}
        s22 = {v.elementLabel: v.data[1] for v in (sf.values if sf else [])}
        s33 = {v.elementLabel: v.data[2] for v in (sf.values if sf else [])}
        # Mises = sqrt(0.5*((s1-s2)^2+(s2-s3)^2+(s3-s1)^2))
        def mises(el):
            a, b, c = s11.get(el, 0.), s22.get(el, 0.), s33.get(el, 0.)
            return ((0.5 * ((a-b)**2 + (b-c)**2 + (c-a)**2)) ** 0.5)

        u3 = {}
        if uf is not None:
            for v in uf.values:
                if hasattr(v, "nodeLabel"):
                    u3[v.nodeLabel] = v.data[2]  # Uz displacement

        S11_hist.append([s11.get(el, float('nan')) for el in el_sorted])
        S22_hist.append([s22.get(el, float('nan')) for el in el_sorted])
        S33_hist.append([s33.get(el, float('nan')) for el in el_sorted])
        SMIS_hist.append([mises(el) for el in el_sorted])
        # average Uz at top z level
        top_z = el_z[el_sorted[-1]]
        top_nodes = [nd.label for nd in inst.nodes if abs(nd.coordinates[2] - top_z) < 1e-4]
        u3_top = sum(u3.get(nd, 0.) for nd in top_nodes) / max(len(top_nodes), 1)
        U3_hist.append(u3_top)

    odb.close()

    return {
        "case": job_name,
        "t":    times,
        "z_mm": [el_z[el] for el in el_sorted],
        "S11":  S11_hist,
        "S22":  S22_hist,
        "S33":  S33_hist,
        "SMISES": SMIS_hist,
        "U3_top": U3_hist,
    }


def main():
    all_data = {}
    for (case_dir_name, job_name) in CASES:
        case_dir = os.path.join(RUN_DIR, case_dir_name)
        print(f"Reading {job_name} ...", file=sys.stderr)
        data = extract(case_dir, job_name)
        all_data[job_name] = data
        print(f"  {len(data['t'])} frames, {len(data['z_mm'])} elements", file=sys.stderr)

    out_path = os.path.expanduser("~/IKM_Hiwi/FEM/visco_2ch_results.json")
    with open(out_path, "w") as f:
        json.dump(to_py(all_data), f, indent=2)
    print(f"Saved: {out_path}", file=sys.stderr)

    for job_name, data in all_data.items():
        t    = data["t"]
        S11  = data["S11"]
        SMIS = data["SMISES"]
        U3   = data["U3_top"]
        z    = data["z_mm"]
        bot, mid, top = 0, len(z)//2, -1
        print(f"\n=== {job_name} ===")
        print(f"  z: bot={z[bot]:.3f}  mid={z[mid]:.3f}  top={z[top]:.3f} mm")
        for i in [0, len(t)//4, len(t)//2, 3*len(t)//4, -1]:
            print(f"  t={t[i]:7.1f}s  S11_bot={S11[i][bot]:+.5f}  "
                  f"S11_top={S11[i][top]:+.5f}  "
                  f"Mises_max={max(SMIS[i]):+.5f}  "
                  f"U3_top={U3[i]:+.5f}  [MPa/mm]")


if __name__ == "__main__":
    main()
