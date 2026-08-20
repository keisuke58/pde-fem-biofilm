#!/bin/sh
# Block any commit whose author/committer identity looks like an AI
# assistant default (e.g. "Claude <noreply@anthropic.com>"), added
# 2026-08-20 after 46 commits in this repo's history turned out to have
# exactly that author -- traced to a Claude Code environment other than
# the primary IKMHIWI03 machine (no global .gitconfig exists there, so it
# wasn't that machine's default -- it came from somewhere else, possibly a
# cloud/remote Claude Code session). See CLAUDE.md and the
# git-history-rewrite-2026-08-20 memory for the full incident.
#
# This only protects commits made from wherever it's installed as
# .git/hooks/pre-commit -- it is NOT installed automatically in a fresh
# clone (git does not track .git/hooks/). Install it in any new clone
# (including after re-cloning the Keio server copy) with:
#
#   cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit
#   chmod +x .git/hooks/pre-commit

check() {
    name="$1"; email="$2"; role="$3"
    case "$email" in
        *anthropic.com*|*noreply@anthropic*)
            echo "pre-commit: refusing -- $role email '$email' looks like an AI-assistant default, not a human identity." >&2
            echo "Set it explicitly: git config user.email \"you@example.com\"" >&2
            exit 1
            ;;
    esac
    case "$name" in
        [Cc]laude|*[Cc]laude\ [Cc]ode*)
            echo "pre-commit: refusing -- $role name '$name' looks like an AI-assistant default, not a human identity." >&2
            echo "Set it explicitly: git config user.name \"Your Name\"" >&2
            exit 1
            ;;
    esac
}

author_name=$(git var GIT_AUTHOR_IDENT | sed 's/ <.*//')
author_email=$(git var GIT_AUTHOR_IDENT | sed -n 's/.*<\(.*\)>.*/\1/p')
committer_name=$(git var GIT_COMMITTER_IDENT | sed 's/ <.*//')
committer_email=$(git var GIT_COMMITTER_IDENT | sed -n 's/.*<\(.*\)>.*/\1/p')

check "$author_name" "$author_email" "author"
check "$committer_name" "$committer_email" "committer"

exit 0
