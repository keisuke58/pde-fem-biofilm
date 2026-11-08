#!/usr/bin/env bash
# bash meetings/keio_deformation/build.sh -> build/deform_ja.pptx (equations rendered from eq/*.tex)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; mkdir -p "$HERE/build/eq"
for f in eq30 flow; do (cd "$HERE/build/eq" && pdflatex -interaction=nonstopmode "$HERE/eq/$f.tex" >/dev/null && pdftoppm -png -r 500 -singlefile "$f.pdf" "$f"); done
NM="$HERE/../2026-11_coupled_examples/build/node_modules"
NODE_PATH="$NM" node "$HERE/deck_deform_ja.js"
