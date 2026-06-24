#!/bin/bash
#PBS -N crown_hconv
#PBS -l nodes=1:ppn=2
#PBS -l walltime=18:00:00
#PBS -q default
#PBS -j oe
#PBS -o /home/nishioka/IKM_Hiwi/FEM/tier2b_real/crown_hconv_${PBS_JOBID}.log
#PBS -m n

# PBS wrapper for the crown h-convergence sweep (run_crown_hconvergence.sh).
# Runs on a COMPUTE node — the bare run script must never be launched on the fifawc login node
# (Abaqus solves are heavy; login-node compute is forbidden on this shared cluster).
#
# 4 LCs x 2 jobs (crown/generic) x 2 orders (C3D4/C3D10) = 16 serial solves (~37 min each, ~10 h).
# walltime 18 h leaves headroom. ppn=2 is courteous on the shared default queue; the solves are
# cpus=1 serial, the extra core only helps gmsh meshing / extraction.
#
# Submit:  qsub /home/nishioka/IKM_Hiwi/FEM/tier2b_real/run_crown_hconvergence_pbs.sh
# Output:  hconv_results.jsonl  (then: python fig_crown_hconvergence.py / crown_iso14801_validation.py)
set -euo pipefail
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
exec bash run_crown_hconvergence.sh
