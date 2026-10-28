#!/usr/bin/env bash
# claude_remote.sh -- start Claude Code Remote Control inside tmux (Keio server, fifa).
#
#   bash scripts/claude_remote.sh          # start, or attach if already running
#   bash scripts/claude_remote.sh stop     # stop the session
#   bash scripts/claude_remote.sh status
#
# Needs a claude.ai login (run `claude` once and /login). The session lives in
# tmux, so closing ssh does not stop it. Detach with Ctrl-b d.
set -u
SESSION=claude-rc
cd "$(git rev-parse --show-toplevel)" || exit 1

command -v tmux >/dev/null || { echo "tmux not found (install it, or run 'claude remote-control' by hand)"; exit 1; }
command -v claude >/dev/null || { echo "claude not found on PATH"; exit 1; }

case "${1:-start}" in
  stop)   tmux kill-session -t "$SESSION" 2>/dev/null && echo "stopped $SESSION" || echo "not running" ;;
  status) tmux has-session -t "$SESSION" 2>/dev/null && echo "running ($SESSION)" || echo "not running" ;;
  start)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "attaching to $SESSION"
    else
      tmux new-session -d -s "$SESSION" -c "$PWD" "claude remote-control"
      echo "started $SESSION"
    fi
    tmux attach -t "$SESSION" ;;
  *) echo "usage: $0 [start|stop|status]"; exit 1 ;;
esac
