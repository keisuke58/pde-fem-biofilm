#!/usr/bin/env bash
# Real-Abaqus single-element check of the PRODUCTION thesis UMAT
# umat_klempt_voigt.f (Option A). Its elastic kernel + analytic Jaumann
# tangent are identical to umat_klempt2025.f (Option C), so this covers
# both. Growth via PREDEF(1)=alpha vs a clamped face, NLGEOM.
# Expect: COMPLETED, maxEquilIter<=2, 0 cutbacks, increment auto-grows
# (= Abaqus confirms the hand-derived consistent tangent is correct).
set -e
cd "$(dirname "$0")"
AB=/home/nishioka/IKM_Hiwi/nife/masterarbeit_ansys_fem/coupling_prototype/abaqus
abaqus job=kv1 input=klempt_voigt_1elem.inp user=$AB/umat_klempt_voigt.f \
    cpus=1 interactive ask_delete=OFF
grep -E "COMPLETED SUCCESSFULLY" kv1.sta && echo "PASS: converged"
