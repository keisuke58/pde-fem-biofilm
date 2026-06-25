#!/usr/bin/env bash
# Real-Abaqus single-element verification of ../../umat_biofilm_visco.f.
# Growth (TEMP=alpha) against a clamped face + viscoelastic relaxation,
# UNSYMM=YES. Confirms the exact consistent tangent gives clean implicit
# convergence (expect maxEquilIter<=2, zero cutbacks, increment auto-grows).
set -e
cd "$(dirname "$0")"
abaqus job=umat_visco_1elem input=umat_visco_1elem.inp \
    user=../../umat_biofilm_visco.f cpus=1 interactive ask_delete=OFF
grep -E "COMPLETED SUCCESSFULLY" umat_visco_1elem.sta && echo "PASS: converged"
