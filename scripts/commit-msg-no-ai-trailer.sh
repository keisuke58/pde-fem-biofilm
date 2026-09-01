#!/bin/sh
# Block any commit MESSAGE that attributes the work to an AI assistant.
#
# This is the companion to pre-commit-no-ai-identity.sh, and it exists
# because that hook has a gap: it checks the author and committer
# identity, but a `Co-Authored-By: Claude <noreply@anthropic.com>`
# trailer is neither -- and GitHub counts co-authors as contributors, so
# such a trailer puts Claude on the repository's Contributors page just
# as surely as a wrong author would. The identity hook cannot see it,
# because pre-commit does not receive the commit message; only
# commit-msg does.
#
# Install (in this clone and in any fresh one):
#
#   cp scripts/commit-msg-no-ai-trailer.sh .git/hooks/commit-msg
#   chmod +x .git/hooks/commit-msg
#
# See CLAUDE.md for the 2026-08-20 incident this line of defence is for.

msg_file="$1"
[ -f "$msg_file" ] || exit 0

# Ignore comment lines -- git's own template mentions nothing relevant,
# but a user's might, and a comment never reaches the stored message.
body=$(grep -v '^#' "$msg_file")

fail() {
    echo "commit-msg: refusing -- $1" >&2
    echo "This repository does not mark commits as AI-assisted; see CLAUDE.md." >&2
    exit 1
}

echo "$body" | grep -qiE '^[[:space:]]*co-authored-by:.*(anthropic\.com|claude)' \
    && fail "the message carries an AI co-author trailer, which would add that identity to the Contributors page."

echo "$body" | grep -qiE '^[[:space:]]*(generated with|created with|written by).*(claude|anthropic)' \
    && fail "the message attributes authorship to an AI assistant."

echo "$body" | grep -qiE '(🤖|:robot:).*(claude|generated)' \
    && fail "the message carries an AI-generation marker."

exit 0
