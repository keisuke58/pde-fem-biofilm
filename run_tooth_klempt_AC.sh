#!/usr/bin/env bash
# run_tooth_klempt_AC.sh
# Mode A が全条件 PASS したら自動で mode C を起動する連鎖ランナー。
# 既に mode A が走っている場合もここから実行すれば多重起動にならない。
#
# Usage:
#   bash run_tooth_klempt_AC.sh           # nohup 推奨
#   nohup bash run_tooth_klempt_AC.sh > run_AC.log 2>&1 &

set -u
cd "$(dirname "$0")"

N_COND=4
RESULTS_A="tooth_klempt_results_A.csv"
POLL=30   # 30秒ごとにチェック

echo "=========================================="
echo "  Klempt AC chain runner — $(date)"
echo "  Waiting for mode A to produce $RESULTS_A with $N_COND PASS/ERROR lines..."
echo "=========================================="

# ── Wait for mode A ─────────────────────────────────────────────────────────
while true; do
    if [ -f "$RESULTS_A" ]; then
        # ヘッダー除く行数 = 完了ジョブ数
        done=$(tail -n +2 "$RESULTS_A" | grep -c "," 2>/dev/null || echo 0)
        pass=$(tail -n +2 "$RESULTS_A" | grep -c ",PASS," 2>/dev/null || echo 0)
        err=$(tail -n +2  "$RESULTS_A" | grep -c ",ERROR," 2>/dev/null || echo 0)
        echo "  [$(date '+%H:%M:%S')]  completed=$done/$N_COND  PASS=$pass  ERROR=$err"

        if [ "$done" -ge "$N_COND" ]; then
            echo ""
            echo "  Mode A finished. Results:"
            cat "$RESULTS_A"
            break
        fi
    else
        echo "  [$(date '+%H:%M:%S')]  $RESULTS_A not yet created..."
    fi
    sleep "$POLL"
done

# ── Launch mode C ────────────────────────────────────────────────────────────
echo ""
echo "=========================================="
echo "  Starting mode C — $(date)"
echo "=========================================="

# inp ファイルが存在しない場合は生成
if [ ! -f "p23_klempt_C_commensal_static.inp" ]; then
    echo "  Generating mode C inp files..."
    python gen_tooth_klempt_umat_inp.py --mode C
fi

bash run_tooth_klempt.sh --mode C

echo ""
echo "=========================================="
echo "  Both modes completed — $(date)"
echo "  Mode A: $RESULTS_A"
echo "  Mode C: tooth_klempt_results_C.csv"
echo ""
echo "  Compare results:"
echo "    echo '=== A ===' && cat tooth_klempt_results_A.csv"
echo "    echo '=== C ===' && cat tooth_klempt_results_C.csv"
echo "=========================================="
