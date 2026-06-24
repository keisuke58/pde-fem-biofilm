#!/bin/bash
# T1 — per-design bone-loss sweep -> crest(L) per geometry -> b_crit(design).
# Reuses periimplantitis_coupon.py (expose x crown_h) + extract_pimp.py (crest_p95), the machinery
# behind run_pimp.sh / run_crown_design.sh. Crown height (occlusal moment arm) is the design knob that
# the peri-implantitis coupon exposes (diameter is NOT parameterised there -> future builder change).
#
# Grid: crown_h {0,4,8} x bone-loss expose {1,3,5,7} = 12 solves (~35 min each, ~7 h).
# Output: t1_design_loss.jsonl  (tag, crown_h, expose, crest_p95, ...). Post: b_crit per crown via the
# vicious-cycle model (rank_designs_bcrit.py, run locally afterwards).
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
CROWNS=(0 4 8)
EXPOSES=(1 3 5 7)

echo "=== build INPs (gmsh_env) ==="
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
rm -f t1_design_loss.jsonl
for ch in "${CROWNS[@]}"; do
  for e in "${EXPOSES[@]}"; do
    tag="t1_h${ch}_e${e}"
    python periimplantitis_coupon.py "$e" bond 2 "$tag" "$ch" 0 2>&1 | grep -E "pimp $tag|wrote|Error" || true
  done
done
conda deactivate 2>/dev/null || true

echo "=== solve + extract (Abaqus) ==="
for ch in "${CROWNS[@]}"; do
  for e in "${EXPOSES[@]}"; do
    tag="t1_h${ch}_e${e}"
    zbt=$(python3 -c "print(10.0-$e)")
    rm -f "pimp_${tag}.lck"
    "$ABQ" job="pimp_${tag}" cpus=1 interactive >"pimp_${tag}.log" 2>&1 || true
    if grep -q COMPLETED "pimp_${tag}.log"; then
      line=$("$ABQ" python extract_pimp.py "$tag" "$zbt" 2>/dev/null | tail -1)
      python3 -c "import json; r=json.loads('''$line'''); r['crown_h']=$ch; r['expose']=$e; open('t1_design_loss.jsonl','a').write(json.dumps(r)+'\n'); print('  h=%smm L=%smm crest_p95=%.1f'%($ch,$e,r['crest_p95']))" || echo "  WARN: parse failed $tag"
    else
      echo "$tag FAILED"; grep -i error "pimp_${tag}.dat" 2>/dev/null | head -2
    fi
  done
done
echo "=== T1 done -> t1_design_loss.jsonl ; next: python rank_designs_bcrit.py ==="
