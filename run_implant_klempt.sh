#!/usr/bin/env bash
# run_implant_klempt.sh — Serial Abaqus runs for helical implant screw Klempt UMAT v3
# 4 conditions, mode A (umat_klempt_voigt.f, phi^2-gated E, PREDEF 1-7)
# Mesh: 48x4x60 C3D8, 3 helical turns, 14640 nodes, 11520 elems
#
# Pre-requisites:
#   python gen_implant_klempt_inp.py --cond {cond} --mode A  (all 4 already generated)
#
# Usage:
#   bash run_implant_klempt.sh           # all 4 conditions
#   bash run_implant_klempt.sh DH_only   # dysbiotic_hobic only (quick test)

set -u
export PATH="$HOME/DassaultSystemes/SIMULIA/Commands:$PATH"
cd "$(dirname "$0")"

UMAT="/home/nishioka/IKM_Hiwi/nife/masterarbeit_ansys_fem/coupling_prototype/abaqus/umat_klempt_voigt.f"
CONDITIONS=(commensal_static commensal_hobic dysbiotic_static dysbiotic_hobic)

if echo "$@" | grep -q "DH_only"; then
    CONDITIONS=(dysbiotic_hobic)
fi

TIMEOUT=3600   # 1 h per job (implant mesh smaller than tooth: 11520 vs ~50k elems)
POLL=10
RESULTS="implant_klempt_results_A.csv"
echo "cond,status,walltime_s" > "$RESULTS"

run_one() {
    local cond="$1"
    local job="p23imp_klempt_A_${cond}"
    local t_start=$(date +%s)

    echo ""
    echo "========================================================"
    echo "  Running: $job"
    echo "  UMAT:    $UMAT"
    echo "  Start:   $(date)"
    echo "========================================================"

    rm -f "${job}".{dat,msg,sta,prt,odb,lck,stt,com} 2>/dev/null
    rm -rf "${job}".simdir 2>/dev/null

    abaqus job="$job" user="$UMAT" cpus=1 ask_delete=OFF >/dev/null 2>&1

    local t=0 status="TIMEOUT"
    while [ $t -lt $TIMEOUT ]; do
        sleep "$POLL"; t=$((t + POLL))
        if grep -q "COMPLETED SUCCESSFULLY" "${job}.sta" 2>/dev/null; then
            status="PASS"; break
        fi
        if grep -qiE "\*\*\*ERROR|HAS NOT BEEN COMPLETED" "${job}.dat" 2>/dev/null; then
            status="ERROR"; break
        fi
        if grep -qiE "\*\*\*ERROR" "${job}.msg" 2>/dev/null; then
            status="ERROR"; break
        fi
        if [ $((t % 60)) -eq 0 ]; then
            echo "  [${t}s] $(tail -1 ${job}.sta 2>/dev/null)"
        fi
    done

    local wall=$(( $(date +%s) - t_start ))
    echo "  → $status  (${wall}s)"
    echo "$cond,$status,$wall" >> "$RESULTS"

    if [ "$status" = "ERROR" ]; then
        grep -A5 "\*\*\*ERROR" "${job}.msg" 2>/dev/null | head -20
        grep -A5 "\*\*\*ERROR" "${job}.dat" 2>/dev/null | head -20
    fi
}

echo "Klempt implant FEM — mode A (phi^2-gated)  $(date)"
echo "UMAT:  $UMAT"
echo "Conds: ${CONDITIONS[*]}"

for cond in "${CONDITIONS[@]}"; do
    run_one "$cond"
done

echo ""
echo "========================================================"
echo "All implant jobs done:"
cat "$RESULTS"
echo "========================================================"
