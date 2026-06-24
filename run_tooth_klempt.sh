#!/usr/bin/env bash
# run_tooth_klempt.sh — Serial Abaqus runs for conformal tooth Klempt UMAT (4 conditions)
#
# IMPORTANT: QSD license has only 50 tokens → strictly serial.
#
# Modes (--mode):
#   A  : umat_klempt_voigt.f   (Option A+B: Voigt E + condition-specific α)
#   C  : umat_klempt2025.f     (Option C: Klempt 2025 multi-species per-species α_s)
#
# Pre-requisites:
#   python JAXFEM/klempt_pde_multispecies.py    (generates klempt_alpha_final_*.npy)
#   python gen_tooth_klempt_umat_inp.py --mode A  (or --mode C)
#
# Usage:
#   bash run_tooth_klempt.sh                    # mode A, all 4 conditions
#   bash run_tooth_klempt.sh --mode C           # mode C, all 4 conditions
#   bash run_tooth_klempt.sh --mode A DH_only   # mode A, dysbiotic_hobic only (test)

set -u
export PATH="$HOME/DassaultSystemes/SIMULIA/Commands:$PATH"
cd "$(dirname "$0")"

UMAT_DIR="/home/nishioka/IKM_Hiwi/nife/masterarbeit_ansys_fem/coupling_prototype/abaqus"

# Parse --mode argument (default A)
MODE="A"
for arg in "$@"; do
    if [ "$arg" = "--mode" ]; then
        MODE_NEXT=1
    elif [ "${MODE_NEXT:-0}" = "1" ]; then
        MODE="$arg"
        MODE_NEXT=0
    fi
done

case "$MODE" in
    A) UMAT="${UMAT_DIR}/umat_klempt_voigt.f"  ;;
    C) UMAT="${UMAT_DIR}/umat_klempt2025.f"    ;;
    *) echo "Unknown mode: $MODE (use A or C)"; exit 1 ;;
esac

TIMEOUT=7200   # 2 h per job (conformal tet mesh is larger than felix square mesh)
POLL=10        # seconds between .sta polls

CONDITIONS=(commensal_static commensal_hobic dysbiotic_static dysbiotic_hobic)

# Optional: run only one condition for a quick test
if echo "$@" | grep -q "DH_only"; then
    CONDITIONS=(dysbiotic_hobic)
fi

results_csv="tooth_klempt_results_${MODE}.csv"
echo "cond,status,walltime_s" > "$results_csv"

run_one() {
    local cond="$1"
    local job="p23_klempt_${MODE}_${cond}"
    local t_start=$(date +%s)

    echo ""
    echo "========================================================"
    echo "  Running: $job"
    echo "  UMAT:    $UMAT"
    echo "  Start:   $(date)"
    echo "========================================================"

    # Clean previous run artifacts
    rm -f "${job}".{dat,msg,sta,prt,odb,lck,stt,com} 2>/dev/null
    rm -rf "${job}".simdir 2>/dev/null

    # Launch detached (no `interactive` — hangs without TTY)
    abaqus job="$job" user="$UMAT" cpus=1 ask_delete=OFF >/dev/null 2>&1

    # Poll .sta for completion or error
    local t=0
    local status="TIMEOUT"
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

        # Progress: print .sta tail every 5 polls
        if [ $((t % 50)) -eq 0 ]; then
            echo "  [${t}s] $(tail -1 ${job}.sta 2>/dev/null)"
        fi
    done

    local t_end=$(date +%s)
    local wall=$((t_end - t_start))
    echo "  → $status  (${wall}s)"
    echo "$cond,$status,$wall" >> "$results_csv"

    if [ "$status" = "ERROR" ]; then
        echo "  --- ERROR excerpt from .msg ---"
        grep -A3 "\*\*\*ERROR" "${job}.msg" 2>/dev/null | head -20
        echo "  --- ERROR excerpt from .dat ---"
        grep -A3 "\*\*\*ERROR" "${job}.dat" 2>/dev/null | head -20
    fi
}

echo "Klempt UMAT tooth FEM — mode=${MODE}  $(date)"
echo "UMAT:       $UMAT"
echo "Conditions: ${CONDITIONS[*]}"

for cond in "${CONDITIONS[@]}"; do
    run_one "$cond"
done

echo ""
echo "========================================================"
echo "All jobs done. Summary:"
cat "$results_csv"
echo ""
echo "Compare results with: python compare_tooth_klempt.py"
echo "========================================================"
