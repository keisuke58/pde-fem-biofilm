#!/usr/bin/env bash
# setup_keio_server.sh -- one-time setup of a clone on the Keio Linux server (fifa).
#
#   bash scripts/setup_keio_server.sh          # from the repository root
#
# Sets the git identity and the two hooks (CLAUDE.md, "Git"), checks Python,
# numpy/scipy and the Abaqus command, and prints what to read next
# (KEIO_SERVER_HANDOFF.ja.md). Safe to run again; it changes only .git/config
# and .git/hooks of this clone.
set -u
cd "$(git rev-parse --show-toplevel)" || exit 1

ok()   { printf '  [ok]   %s\n' "$1"; }
warn() { printf '  [warn] %s\n' "$1"; }

echo "git"
git config --local user.name  "keisuke nishioka"
git config --local user.email "128669518+keisuke58@users.noreply.github.com"
git config --local commit.gpgsign false
ok "identity: $(git var GIT_AUTHOR_IDENT | sed 's/ [0-9].*//')"
cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit
cp scripts/commit-msg-no-ai-trailer.sh  .git/hooks/commit-msg  && chmod +x .git/hooks/commit-msg
ok "hooks pre-commit and commit-msg installed"
if git log -1 --format=%an 2>/dev/null | grep -qi claude; then
  warn "last commit author looks like an AI identity; see CLAUDE.md before committing"
fi
branch=$(git rev-parse --abbrev-ref HEAD)
ok "branch: $branch"

echo "python"
if command -v python3 >/dev/null; then
  ok "$(python3 --version 2>&1)"
  python3 -c "import numpy, scipy; print('  [ok]   numpy', numpy.__version__, 'scipy', scipy.__version__)" \
    || warn "numpy/scipy missing: python3 -m pip install --user numpy scipy matplotlib"
else
  warn "python3 not found"
fi

echo "abaqus"
abq=${ABAQUS:-abaqus}
if command -v "$abq" >/dev/null; then
  ok "$abq found at $(command -v "$abq")"
  echo "         next: $abq verify -user_std   (user subroutines compile?)"
else
  warn "'$abq' not on PATH; set ABAQUS=/path/to/abaqus (run_comp.sh reads it)"
fi
echo "         work dir: ${WORKROOT:-$HOME/abaqus_work} (set WORKROOT for a larger disk)"

echo
echo "Next: read KEIO_SERVER_HANDOFF.ja.md (sections 3-6), then KEIO_PLAN.ja.md section 0."
echo "First run: section 4 of the hand-off (wp2_imp_free_ph01, compare with results_1007/)."
