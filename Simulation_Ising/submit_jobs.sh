#!/bin/bash
#SBATCH -N 1
#SBATCH -t 10:00:00

module purge
module load 2020
module load R/4.0.2-intel-2020a

cp -r "$HOME"/GDSimIsing "$TMPDIR"
cd "$TMPDIR"/GDSimIsing

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla simulation.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/GDSimIsing/output