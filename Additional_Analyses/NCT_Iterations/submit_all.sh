#!/bin/bash

mkdir "$TMPDIR"/GDSim_NCT/

cd "$HOME"/GDSim_NCT

sbatch -a 1-50 submit_jobs.sh




