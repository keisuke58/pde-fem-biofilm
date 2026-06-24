#!/bin/bash
#PBS -N t1_design
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -q default
#PBS -j oe
#PBS -o /home/nishioka/IKM_Hiwi/FEM/tier2b_real/t1_design_${PBS_JOBID}.log
#PBS -m n

# PBS wrapper for the T1 per-design bone-loss sweep. Submit chained after the sensitivity job so the
# three heavy Abaqus sweeps run sequentially (shared-cluster courtesy):
#   qsub -W depend=afterany:<sens_jobid> run_t1_design_loss_pbs.sh
cd /home/nishioka/IKM_Hiwi/FEM/tier2b_real
exec bash run_t1_design_loss.sh
