#!/bin/bash

mkdir "$TMPDIR"/GDSimIsing/

cd "$HOME"/GDSimIsing

sbatch -a 1-200 submit_jobs.sh




