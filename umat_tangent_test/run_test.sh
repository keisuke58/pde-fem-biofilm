#!/usr/bin/env bash
# Standalone verification of the consistent tangent in ../umat_biofilm_visco.f
# Compiles the production UMAT against a stub ABA_PARAM.INC (so it builds
# outside Abaqus) plus the driver, then checks DDSDDE vs a central-difference
# reference. Mirrors ../phase2_patch_test.py.
#   PASS: symmetry < 5e-3 and max rel err vs FD < 1e-2 for all cases.
set -e
cd "$(dirname "$0")"
gfortran -O2 -ffixed-line-length-72 -I. \
    ../umat_biofilm_visco.f test_umat_tangent.f -o test_umat
./test_umat
