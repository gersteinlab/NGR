#!/bin/bash

#SBATCH --partition=pi_gerstein,general
#SBATCH --job-name=generate_all_mobility_lists
#SBATCH --mem=30GB
#SBATCH --time=23:59:59

# run with mobility lists section uncommented
python evaluate.py

echo "Done."
date