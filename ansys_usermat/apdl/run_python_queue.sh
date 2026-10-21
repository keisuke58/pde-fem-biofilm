#!/usr/bin/env bash
# run_python_queue.sh -- runs coupled_living_consumption.py jobs one after another and commits + pushes each
# result JSON as soon as it exists (a cloud container can be reclaimed; finished runs survive in git).
#   bash ansys_usermat/apdl/run_python_queue.sh LANE_NAME "args of run 1" "args of run 2" ...
# Base options (the agreed set-up, 21 Oct) are added to every run.
cd "$(dirname "$0")/../.." || exit 1
lane=$1; shift
BASE="--v2 --phi-min 0.1 --front 0.1 --fortran --ansys-bc --start-position"
log=logs_plaque/pyqueue_${lane}.log
mkdir -p logs_plaque
for args in "$@"; do
  echo "[$(date -u +%FT%TZ)] start: $args" >> "$log"
  before=$(ls ansys_usermat/apdl/results/2026-10-3d/*.json 2>/dev/null | sort)
  python ansys_usermat/apdl/coupled_living_consumption.py $BASE $args >> "$log" 2>&1
  rc=$?
  after=$(ls -t ansys_usermat/apdl/results/2026-10-3d/*.json | head -n 3)
  new=$(for f in $after; do [ "$f" -nt "$log.stamp" ] 2>/dev/null && echo "$f"; done)
  echo "[$(date -u +%FT%TZ)] end ($rc): $args" >> "$log"
  for f in $(git status --porcelain ansys_usermat/apdl/results/2026-10-3d/ | awk '{print $2}'); do
    git add "$f"
  done
  if ! git diff --cached --quiet; then
    git commit -q -m "Python queue $lane: $args" && \
    for i in 1 2 3 4; do git pull -q --rebase origin claude/plan-next-hxjjve && git push -q origin HEAD && break; sleep $((2**i)); done
  fi
  touch "$log.stamp"
done
echo "[$(date -u +%FT%TZ)] lane $lane done" >> "$log"
