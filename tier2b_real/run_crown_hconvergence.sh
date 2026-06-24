#!/usr/bin/env bash
# h-convergence sweep on the crowned headline assembly (tier2b_crown / tier2b_generic).
#
# Varies the parametric crown tet edge LC in mesh_crown.py via the CROWN_LC env var (introduced as a
# tiny env-overridable patch — default 0.40 preserves headline behaviour). For each LC the crown mesh
# is re-built, the assembly re-assembled at the SAME headline material/geometry, solved at C3D4
# linear tets, then ALSO converted to C3D10 and re-solved so each LC delivers two points (linear and
# quadratic). The thesis-metric crestal peri-implant bone p95 is extracted into hconv_results.jsonl.
#
# Purpose: complement the C3D10 element-order audit (notes/fem_element_fidelity_2026-06-18.md, §2)
# with the matching h-convergence story so the headline (~275 k C3D4 elements, LC=0.40 crown) can be
# certified resolution-converged at the thesis metric, and the C3D10 jump (+8.9 %) can be separated
# from any residual h-discretisation error via Richardson extrapolation (see fig_crown_hconvergence.py).
#
# Compute budget: each solve ~37 min, 4 LCs × 2 jobs (crown / generic) × 2 orders (linear+quadratic)
# = 16 solves ~10 h serial (QSD=50 forces strict serial). Run on the cluster.
#
#   cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
#   bash run_crown_hconvergence.sh
set -euo pipefail
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real

ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
LCS=(0.55 0.40 0.30 0.22)            # halving-doubling sequence; 0.40 = current headline

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

rm -f hconv_results.jsonl

build_and_solve() {
    local lc=$1
    local tag_lc="lc${lc//./p}"
    local crown_job="tier2b_crown_${tag_lc}"
    local generic_job="tier2b_generic_${tag_lc}"

    echo "=== [LC=$lc] regenerate crown mesh (via mesh_crown_lc.py wrapper) ==="
    python mesh_crown_lc.py "$lc" | tail -2

    for job_pair in "cache_implant.npz $crown_job crown" "cache_implant_generic.npz $generic_job generic"; do
        read -r cache jobname kind <<<"$job_pair"
        echo "--- [$kind LC=$lc] assemble $jobname (C3D4) ---"
        python build_assembly.py "$cache" "$jobname" 2>&1 | tail -2

        echo "--- [$kind LC=$lc] convert to C3D10 ($jobname -> ${jobname}_q) ---"
        python convert_c3d4_to_c3d10.py "${jobname}.inp" "${jobname}_q.inp" 2>&1 | tail -2

        for order_tag in "" "_q"; do
            local J="${jobname}${order_tag}"
            local order=$([ -z "$order_tag" ] && echo C3D4 || echo C3D10)
            rm -f "${J}.lck"
            echo ">>> solve $J ($order)"
            "$ABQ" job="$J" cpus=1 interactive >"${J}.log" 2>&1 || true
            if grep -q COMPLETED "${J}.log"; then
                # extract_tier2b_q.py emits ${J}_field.json with x,y,z,mat,vmg,vmo for BOTH C3D4 and
                # C3D10 inps (extract_pimp_field.py was the wrong extractor — pimp_ prefix, "BONE"
                # elset, r/vm schema — and silently failed, killing the sweep under set -e).
                "$ABQ" python extract_tier2b_q.py "$J" >"${J}.extract.log" 2>&1 || \
                    echo "  WARN: extract failed for $J (see ${J}.extract.log)"
                if [ -f "${J}_field.json" ]; then
                    # one-line p95 via the EXACT thesis metric; never fatal to the sweep (|| true)
                    python3 - "$J" "$lc" "$kind" "$order" <<'PY' || echo "  WARN: p95 failed for $J"
import json, sys, numpy as np
job, lc, kind, order = sys.argv[1:]
T23 = np.array([-69.4, -41.0]); CREST = 29.0
d = json.load(open(f"{job}_field.json"))["els"]
x = np.array([e["x"] for e in d]); y = np.array([e["y"] for e in d]); z = np.array([e["z"] for e in d])
vm = np.array([e["vmo"] for e in d]); mat = np.array([e["mat"] for e in d])
r = np.hypot(x - T23[0], y - T23[1])
m = np.isin(mat, ["CORTICAL", "CANCELLOUS", "BONE"]) & (r <= 3.0) & (z >= CREST - 3) & (z <= CREST + 1.5)
p95 = float(np.percentile(vm[m], 95))
n_el = json.load(open(f"{job}_field.json")).get("n_el", len(d))
rec = {"job": job, "lc": float(lc), "kind": kind, "order": order,
       "n_el": int(n_el), "n_disk_el": int(m.sum()), "crestal_p95_MPa": p95}
print(json.dumps(rec))
open("hconv_results.jsonl", "a").write(json.dumps(rec) + "\n")
PY
                else
                    echo "  WARN: ${J}_field.json missing — skipping p95 for $J (sweep continues)"
                fi
            else
                echo "$J FAILED"
                grep -iE "error|abort" "${J}.dat" 2>/dev/null | head -3 || true
            fi
        done
    done
}

for lc in "${LCS[@]}"; do
    build_and_solve "$lc"
done

conda deactivate 2>/dev/null || true
echo "=== h-convergence sweep done -> hconv_results.jsonl ==="
echo "Next: python fig_crown_hconvergence.py"
