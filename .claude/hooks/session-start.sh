#!/bin/bash
# SessionStart hook — restores what a cloud container does not carry over.
#
# Three problems this exists for, all observed on 2026-09-01:
#
#  1. The container ships a global /root/.gitconfig with user.name=Claude and
#     user.email=noreply@anthropic.com. Any clone without a local user.* picks
#     it up, which is exactly how 46 commits ended up mis-attributed on the
#     Contributors page in August (see CLAUDE.md). A mid-session container
#     restart wiped this repo's local identity and git resolved straight back
#     to the AI default -- the next commit would have repeated the incident.
#
#  2. The anti-attribution hooks live in scripts/ because git does not track
#     .git/hooks/, so a fresh container has neither installed.
#
#  3. A restart loses pip packages and gfortran, so the Fortran-dependent
#     tests silently skip rather than fail -- the worst outcome, since the
#     suite still reports green.
#
# Idempotent and non-interactive; safe to re-run.
set -uo pipefail

cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" || exit 0

log() { echo "[session-start] $*"; }

# --- 1. git identity ---------------------------------------------------
# Set unconditionally rather than only when unset: the failure mode is a
# *wrong* value inherited from the global config, not a missing one.
git config --local user.name  "keisuke nishioka"
git config --local user.email "128669518+keisuke58@users.noreply.github.com"
git config --local commit.gpgsign false   # global points at a key not present here
log "git identity: $(git var GIT_AUTHOR_IDENT | sed 's/ [0-9].*//')"

case "$(git var GIT_AUTHOR_IDENT)" in
    *anthropic.com*|Claude\ *)
        log "ERROR: identity is still an AI default -- do not commit until fixed" ;;
esac

# --- 2. the two guards -------------------------------------------------
for pair in "pre-commit-no-ai-identity.sh:pre-commit" "commit-msg-no-ai-trailer.sh:commit-msg"; do
    src="scripts/${pair%%:*}"; dst=".git/hooks/${pair##*:}"
    if [ -f "$src" ]; then
        cp "$src" "$dst" && chmod +x "$dst" && log "installed $dst"
    fi
done

# --- 3. toolchain ------------------------------------------------------
python -m pip install -q -r requirements.txt 2>&1 | tail -1
# JAX is deliberately absent from requirements.txt; the tangent and ecology
# tests importorskip without it, which means they pass by not running.
python -c "import jax" 2>/dev/null || python -m pip install -q "jax[cpu]" 2>&1 | tail -1

# gfortran gates the UMAT/USERMAT cross-check and every Fortran-driven test.
if ! command -v gfortran >/dev/null 2>&1; then
    log "installing gfortran"
    (sudo apt-get update -qq && sudo apt-get install -y -qq gfortran) >/dev/null 2>&1 \
        || apt-get install -y -qq gfortran >/dev/null 2>&1 \
        || log "WARNING: gfortran unavailable -- Fortran tests will SKIP, not fail"
fi

log "gfortran: $(command -v gfortran || echo MISSING) | jax: $(python -c 'import jax;print(jax.__version__)' 2>/dev/null || echo MISSING)"
log "ready"
