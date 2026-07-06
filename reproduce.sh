#!/usr/bin/env bash
#
# Regenerate the committed figures and validation metrics from source, then run
# the test suite. Requires the deps in requirements.txt (pip install -r ...).
# JAX is not needed for these steps.
#
#   ./reproduce.sh
#
set -euo pipefail
cd "$(dirname "$0")"

run() { echo; echo "==> $*"; "$@"; }

run python3 plot_heine_composition.py          # -> assets/heine_species_composition.png
run python3 plot_heine_phi_psi.py              # -> assets/heine_phi_psi_joint.png
run python3 validate_composition.py            # -> assets/validation_composition_dysbiotic.png + metrics

echo; echo "==> repo audit (runnable subset; non-fatal)"
python3 JAXFEM/audit_all.py --quick || true

run python3 -m pytest

echo; echo "All reproduction steps completed. See assets/ and _validation/."
