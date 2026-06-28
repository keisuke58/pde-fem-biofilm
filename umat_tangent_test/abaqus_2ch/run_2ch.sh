#!/usr/bin/env bash
# Real-Abaqus 1-element verification of umat_biofilm_visco_2ch.f.
# Two-channel Prony (tau1=5s, tau2=50s) + growth (alpha=0.1) vs clamped face.
# UNSYMM=YES.  Expect maxEquilIter<=2, zero cutbacks, auto-growing increments.
set -e
cd "$(dirname "$0")"
abaqus job=umat_visco_2ch_1elem input=umat_visco_2ch_1elem.inp \n    user=../../umat_biofilm_visco_2ch.f cpus=1 interactive ask_delete=OFF
grep -E "COMPLETED SUCCESSFULLY" umat_visco_2ch_1elem.sta && echo "PASS: 2ch UMAT converged"
