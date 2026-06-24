#!/usr/bin/env bash
# Implant-coupon design sweep: build all INPs (gmsh_env), then solve+extract all (Abaqus). Sequential.
# Columns: tag D L pitch taper order angle  (no set -e: a failed solve must not abort the whole sweep)
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus

# tag      D    L   pitch taper order angle
COUPONS=(
 "base   4.1  10  0.8  0.0  2  30"
 "D35    3.5  10  0.8  0.0  2  30"
 "D48    4.8  10  0.8  0.0  2  30"
 "L8     4.1   8  0.8  0.0  2  30"
 "L12    4.1  12  0.8  0.0  2  30"
 "P06    4.1  10  0.6  0.0  2  30"
 "P10    4.1  10  1.0  0.0  2  30"
 "TP3    4.1  10  0.8  0.3  2  30"
 "A0     4.1  10  0.8  0.0  2   0"
 "C3D4   4.1  10  0.8  0.0  1  30"
)

echo "=== Phase 1: build INPs (gmsh_env) ==="
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
rm -f coupon_params.jsonl coupon_results.jsonl
for row in "${COUPONS[@]}"; do
  read tag D L p t o a <<< "$row"
  python implant_coupon.py "$D" "$L" "$p" "$t" "$o" "$a" "$tag" 2>&1 | grep -E "coupon $tag|wrote|Error" || true
done
conda deactivate 2>/dev/null || true

echo "=== Phase 2: solve + extract (Abaqus) ==="
for row in "${COUPONS[@]}"; do
  read tag D L p t o a <<< "$row"
  rm -f coupon_${tag}.lck
  $ABQ job=coupon_${tag} cpus=1 interactive > coupon_${tag}.log 2>&1 || true
  if grep -q "COMPLETED" coupon_${tag}.log; then
    $ABQ python extract_coupon.py "$tag" "$L" 2>/dev/null | tail -1
  else
    echo "$tag FAILED"; grep -i "error" coupon_${tag}.dat 2>/dev/null | head -2
  fi
done
echo "=== sweep done ==="
