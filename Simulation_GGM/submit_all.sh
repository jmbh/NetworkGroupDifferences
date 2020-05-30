#!/bin/bash

mkdir "$TMPDIR"/GDSim/

cd "$HOME"/GDSim

sbatch -a 1-200 submit_jobs.sh