#!/bin/bash
#SBATCH -N 1
#SBATCH -t 04:00:00

module load 2019 Anaconda3
source activate my_root

#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/GDSim "$TMPDIR"
cd "$TMPDIR"/GDSim

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla simulation.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/GDSim/output