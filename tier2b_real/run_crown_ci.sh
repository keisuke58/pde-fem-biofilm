#!/usr/bin/env bash
# Crown-to-implant (C/I) ratio x marginal-bone-loss vicious cycle.
# Sweep marginal bone loss (expose) for a bare abutment (crown_h=0) vs a crowned implant (crown_h=7mm,
# load offset to the occlusal table).  As bone resorbs the lever arm grows AND the anchored length
# shrinks -> the crowned case has a STEEPER crestal-stress-vs-loss feedback A(loss).  Bonded coupon.
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
EXPOSES=(0 1 2 3 4)
CROWNS=(0 7)              # 0 = bare abutment, 7 = +7 mm crown occlusal lever arm

echo "=== build (gmsh_env) ==="
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
rm -f ci_params.jsonl ci_results.jsonl
for ch in "${CROWNS[@]}"; do
  for e in "${EXPOSES[@]}"; do
    tag="ci_c${ch}_e${e}"
    python periimplantitis_coupon.py "$e" bond 2 "$tag" "$ch" 2>&1 | grep -E "pimp $tag|wrote|Error" || true
  done
done
conda deactivate 2>/dev/null || true

echo "=== solve + extract (Abaqus) ==="
for ch in "${CROWNS[@]}"; do
  for e in "${EXPOSES[@]}"; do
    tag="ci_c${ch}_e${e}"
    zbt=$(python3 -c "print(10.0-$e)")
    rm -f pimp_${tag}.lck
    $ABQ job=pimp_${tag} cpus=1 interactive > pimp_${tag}.log 2>&1 || true
    if grep -q "COMPLETED" pimp_${tag}.log; then
      line=$($ABQ python extract_pimp.py "$tag" "$zbt" 2>/dev/null | tail -1)
      # tag back the crown_h + expose for the figure
      python3 -c "import json,sys; r=json.loads('''$line'''); r['crown_h']=$ch; r['expose']=$e; open('ci_results.jsonl','a').write(json.dumps(r)+'\n'); print('  %-12s loss=%smm crown=%smm crest_p95=%.1f MPa'%(r['tag'],$e,$ch,r['crest_p95']))"
    else
      echo "$tag FAILED"; grep -i "error" pimp_${tag}.dat 2>/dev/null | head -2
    fi
  done
done
echo "=== crown-CI sweep done -> ci_results.jsonl ==="
