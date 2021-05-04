#!/bin/bash

mkdir "$TMPDIR"/GDSim_sparsityG/

cd "$HOME"/GDSim_sparsityG

sbatch -a 11-110 submit_jobs.sh


