#!/bin/bash

mkdir "$TMPDIR"/GDSim/

cd "$HOME"/GDSim

sbatch -a 201 submit_jobs.sh