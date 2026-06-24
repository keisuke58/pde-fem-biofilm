#!/usr/bin/env bash
# Material-parameter sensitivity sweep (one-at-a-time, OAT) on the crowned headline assembly.
#
# Targets the soft / uncertain materials in MATS (build_assembly.py): GINGIVA (E ~ 3 MPa with literature
# range 1.5-30 MPa), PDL (~50 MPa, lit 10-100 MPa), CEMENTUM (~15 GPa, lit 8-25 GPa), and BIOFILM
# (~1 MPa for the dysbiotic growth collar). Each is perturbed at 0.5x, 1.0x (baseline echo),
# 2.0x of the headline value, keeping all other materials at headline values.
#
# Output: per-perturbation crestal-p95 (thesis metric) → sens_results.jsonl. fig_crown_sensitivity.py
# turns this into a tornado plot ranking the materials by their dp95/dE.
#
# Compute: 4 mats × 2 non-baseline levels × 2 jobs (crown / generic) = 16 solves, plus 1 baseline
# (already in hconv_results.jsonl if h-convergence has been run). ~10 h serial.
#
#   cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
#   bash run_crown_sensitivity.sh
set -euo pipefail
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real

ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus

# Headline E [MPa] and nu (from build_assembly.py MATS); perturb E only, nu held.
# name        E_baseline  nu     low (E×0.5)   high (E×2.0)
declare -A E_BASE=( [GINGIVA]=3.0  [PDL]=50.0   [CEMENTUM]=15000.0  [BIOFILM]=1.0 )
declare -A NU=(     [GINGIVA]=0.45 [PDL]=0.45   [CEMENTUM]=0.31     [BIOFILM]=0.45 )
LEVELS=(0.5 2.0)   # 1.0 = baseline (use the standard headline result instead of re-solving)

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

rm -f sens_results.jsonl

extract_one() {
    local job=$1 mat=$2 mul=$3 kind=$4
    # extract_tier2b_q.py emits {job}_field.json (x,y,z,mat,vmg,vmo) for C3D4 & C3D10; the old
    # extract_pimp_field.py (pimp_ prefix / "BONE" elset / r,vm schema) silently produced nothing.
    "$ABQ" python extract_tier2b_q.py "$job" >"${job}.extract.log" 2>&1 || \
        echo "  WARN: extract failed for $job (see ${job}.extract.log)"
    [ -f "${job}_field.json" ] || { echo "  WARN: ${job}_field.json missing — skip $job"; return 0; }
    python3 - "$job" "$mat" "$mul" "$kind" <<'PY' || echo "  WARN: p95 failed for $job"
import json, sys, numpy as np
job, mat, mul, kind = sys.argv[1:]
T23 = np.array([-69.4, -41.0]); CREST = 29.0
d = json.load(open(f"{job}_field.json"))["els"]
x = np.array([e["x"] for e in d]); y = np.array([e["y"] for e in d]); z = np.array([e["z"] for e in d])
vm = np.array([e["vmo"] for e in d]); m = np.array([e["mat"] for e in d])
r = np.hypot(x - T23[0], y - T23[1])
sel = np.isin(m, ["CORTICAL", "CANCELLOUS", "BONE"]) & (r <= 3.0) & (z >= CREST - 3) & (z <= CREST + 1.5)
p95 = float(np.percentile(vm[sel], 95))
rec = {"job": job, "mat": mat, "mul": float(mul), "kind": kind, "crestal_p95_MPa": p95}
print(json.dumps(rec))
open("sens_results.jsonl", "a").write(json.dumps(rec) + "\n")
PY
}

for mat in "${!E_BASE[@]}"; do
    E0=${E_BASE[$mat]}; nu=${NU[$mat]}
    for mul in "${LEVELS[@]}"; do
        E=$(python3 -c "print($E0 * $mul)")
        tag="sens_${mat,,}_m${mul//./p}"
        override=$(python3 -c "import json; print(json.dumps({'$mat': [$E, $nu]}))")

        echo "=== [$mat × $mul → E=$E MPa] ==="
        for spec in "cache_implant.npz tier2b_crown_${tag} crown" \
                    "cache_implant_generic.npz tier2b_generic_${tag} generic"; do
            read -r cache jobname kind <<<"$spec"
            MATS_OVERRIDE="$override" python build_assembly_override.py "$cache" "$jobname" 2>&1 | tail -2
            rm -f "${jobname}.lck"
            "$ABQ" job="$jobname" cpus=1 interactive >"${jobname}.log" 2>&1 || true
            if grep -q COMPLETED "${jobname}.log"; then
                extract_one "$jobname" "$mat" "$mul" "$kind"
            else
                echo "$jobname FAILED"
            fi
        done
    done
done

conda deactivate 2>/dev/null || true
echo "=== material-sensitivity sweep done -> sens_results.jsonl ==="
echo "Next: python fig_crown_sensitivity.py"
