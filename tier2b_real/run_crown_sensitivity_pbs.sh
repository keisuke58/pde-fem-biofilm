#!/bin/bash
#PBS -N crown_sens
#PBS -l nodes=1:ppn=2
#PBS -l walltime=18:00:00
#PBS -q default
#PBS -j oe
#PBS -o /home/nishioka/IKM_Hiwi/FEM/tier2b_real/crown_sens_${PBS_JOBID}.log
#PBS -m n

# PBS wrapper for the OAT material-sensitivity sweep (run_crown_sensitivity.sh).
# Runs on a COMPUTE node. Submit with a dependency on the h-convergence job so the two heavy Abaqus
# sweeps never run concurrently on the shared cluster (courtesy):
#   qsub -W depend=afterany:<hconv_jobid> run_crown_sensitivity_pbs.sh
#
# 4 soft materials (GINGIVA, PDL, CEMENTUM, BIOFILM) x 2 levels (0.5x, 2.0x) x 2 jobs (crown/generic)
# = 16 solves; uses extract_tier2b_q.py (extractor fix). Output: sens_results.jsonl
# Next: python fig_crown_sensitivity.py  -> tornado plot of dp95/dE per material.
set -euo pipefail
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
exec bash run_crown_sensitivity.sh
