#!/usr/bin/env bash
# Study (2): restoration DESIGN sweep -- crown height (moment arm) is the design knob.
# Bonded ISO-14801 coupon, fixed bone level (expose=1 mm), crown height varied 0..10 mm.
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
ABQ=~/DassaultSystemes/SIMULIA/Commands/abaqus
HEIGHTS=(0 2 4 6 8 10); EXPOSE=1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmsh_env; export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
rm -f design_results.jsonl
for h in "${HEIGHTS[@]}"; do
  tag="dz_h${h}"
  python periimplantitis_coupon.py "$EXPOSE" bond 2 "$tag" "$h" 0 2>&1 | grep -E "pimp $tag" || true
done
conda deactivate 2>/dev/null || true
zbt=$(python3 -c "print(10.0-$EXPOSE)")
for h in "${HEIGHTS[@]}"; do
  tag="dz_h${h}"; rm -f pimp_${tag}.lck
  $ABQ job=pimp_${tag} cpus=1 interactive > pimp_${tag}.log 2>&1 || true
  if grep -q COMPLETED pimp_${tag}.log; then
    line=$($ABQ python extract_pimp.py "$tag" "$zbt" 2>/dev/null | tail -1)
    python3 -c "import json; r=json.loads('''$line'''); r['crown_h']=$h; open('design_results.jsonl','a').write(json.dumps(r)+'\n'); print('  h=%smm crest_p95=%.1f ti_max=%.1f disp=%.0fum'%($h,r['crest_p95'],r['ti_max'],r['disp_um']))"
  else echo "$tag FAILED"; fi
done
echo "=== crown-height design sweep done ==="
