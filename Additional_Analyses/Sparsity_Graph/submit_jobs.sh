#!/bin/bash
#SBATCH -N 1
#SBATCH -t 02:00:00

module purge
module load 2020
module load R/4.0.2-intel-2020a

cp -r "$HOME"/GDSim_sparsityG "$TMPDIR"
cd "$TMPDIR"/GDSim_sparsityG

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla sim_sparsityG.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/GDSim_sparsityG/output