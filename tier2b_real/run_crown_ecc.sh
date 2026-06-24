#!/usr/bin/env bash
# Study (3): ECCENTRIC occlusal contact -> saucerisation direction.
# Crowned coupon (crown_h=7), fixed bone level, lateral bite offset (+x = buccal) varied 0..4 mm.
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
ECCS=(0 1 2 3 4); EXPOSE=1; CROWN_H=7
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env; export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
rm -f ci_ecc_results.jsonl
for ec in "${ECCS[@]}"; do
  tag="ecc_${ec}"
  python periimplantitis_coupon.py "$EXPOSE" bond 2 "$tag" "$CROWN_H" "$ec" 2>&1 | grep -E "pimp $tag" || true
done
conda deactivate 2>/dev/null || true
zbt=$(python3 -c "print(10.0-$EXPOSE)")
for ec in "${ECCS[@]}"; do
  tag="ecc_${ec}"; rm -f pimp_${tag}.lck
  $ABQ job=pimp_${tag} cpus=1 interactive > pimp_${tag}.log 2>&1 || true
  if grep -q COMPLETED pimp_${tag}.log; then
    $ABQ python extract_pimp_angular.py "$tag" "$zbt" 12 2>/dev/null | tail -1
  else echo "$tag FAILED"; fi
done
echo "=== eccentric saucerisation sweep done ==="
