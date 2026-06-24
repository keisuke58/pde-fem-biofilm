#!/usr/bin/env bash
# Study (1): progressive marginal bone loss. Bonded coupon, bone level lowered in steps.
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
# tag expose(mm) iface order  (bone_top z = 10 - expose)
ROWS=("e2 2 bond 2" "e4 4 bond 2" "e6 6 bond 2" "e8 8 bond 2")

echo "=== build (gmsh_env) ==="
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
rm -f pimp_params.jsonl pimp_results.jsonl
for r in "${ROWS[@]}"; do
  read tag e iface o <<< "$r"
  python periimplantitis_coupon.py "$e" "$iface" "$o" "$tag" 2>&1 | grep -E "pimp $tag|wrote|Error" || true
done
conda deactivate 2>/dev/null || true

echo "=== solve + extract (Abaqus) ==="
for r in "${ROWS[@]}"; do
  read tag e iface o <<< "$r"
  zbt=$(python3 -c "print(10.0-$e)")
  rm -f pimp_${tag}.lck
  $ABQ job=pimp_${tag} cpus=1 interactive > pimp_${tag}.log 2>&1 || true
  if grep -q "COMPLETED" pimp_${tag}.log; then
    $ABQ python extract_pimp.py "$tag" "$zbt" 2>/dev/null | tail -1
  else
    echo "$tag FAILED"; grep -i "error" pimp_${tag}.dat 2>/dev/null | head -2
  fi
done
echo "=== pimp progression done ==="
