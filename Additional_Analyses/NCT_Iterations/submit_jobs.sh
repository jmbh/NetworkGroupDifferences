#!/bin/bash
#SBATCH -N 1
#SBATCH -t 05:00:00

module purge
module load 2020
module load R/4.0.2-intel-2020a

cp -r "$HOME"/GDSim_NCT "$TMPDIR"
cd "$TMPDIR"/GDSim_NCT

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla NCT_sim.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/GDSim_NCT/output