#!/usr/bin/env bash
# Build the decks of the two coupled examples (separate and mixed start, with antibiotic).
#   bash meetings/2026-11_coupled_examples/build.sh [en|ja|all]
# Needs: pdflatex (mathptmx), pdftoppm, node with pptxgenjs (npm install pptxgenjs in build/ if missing);
# the videos need pyvista + libosmesa6 + imageio-ffmpeg (ansys_usermat/apdl/anim_coupled_tvd_3d.py).
# Output: build/coupled_examples_en.pptx, build/report_{separate,mixed}_ab_academic.pptx (not in git).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"; B="$HERE/build"
mkdir -p "$B/eq"
for f in eqs point vm gci pe; do
  (cd "$B/eq" && pdflatex -interaction=nonstopmode "$HERE/eq/$f.tex" >/dev/null && pdftoppm -png -r 600 -singlefile "$f.pdf" "$f")
done
pdftoppm -png -r 600 -singlefile "$ROOT/assets/fig_coupled_3d_panels_separate_ab_alpha.pdf" "$B/s12_coupled_3d_separate_ab"
pdftoppm -png -r 600 -singlefile "$ROOT/assets/fig_coupled_3d_panels_mixed_ab_alpha.pdf" "$B/s12m_coupled_3d_mixed_ab"
pdftoppm -png -r 400 -singlefile "$ROOT/assets/fig_obstacle_columns_3d.pdf" "$B/geom_obst3d"
pdftoppm -png -r 400 -singlefile "$ROOT/assets/fig_maze_k24_h4.pdf" "$B/geom_maze"
for c in separate mixed; do
  [ -f "$ROOT/assets/fig_setup_${c}_ab.png" ] || python3 "$ROOT/ansys_usermat/apdl/fig_setup_separate_ab.py" $([ $c = mixed ] && echo mixed)
  [ -f "$ROOT/assets/anim_coupled_tvd_${c}_ab_inhliv_transc_5panel.mp4" ] || python3 "$ROOT/ansys_usermat/apdl/anim_coupled_tvd_3d.py" ${c}_ab inhliv transc
done
[ -d "$B/node_modules/pptxgenjs" ] || (cd "$B" && npm install --silent pptxgenjs >/dev/null)
export NODE_PATH="$B/node_modules"
case "${1:-all}" in
  en) node "$HERE/deck_en.js" ;;
  ja) node "$HERE/deck_ja.js" separate; node "$HERE/deck_ja.js" mixed ;;
  geom) node "$HERE/deck_geom_en.js" ;;
  k26c4) node "$HERE/deck_k26c4_en.js" ;;
  *)  node "$HERE/deck_en.js"; node "$HERE/deck_ja.js" separate; node "$HERE/deck_ja.js" mixed; node "$HERE/deck_geom_en.js"; node "$HERE/deck_k26c4_en.js" ;;
esac
