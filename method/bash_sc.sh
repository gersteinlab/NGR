#!/bin/bash

#SBATCH --time=47:59:00
#SBATCH --mem=30Gb
#SBATCH --partition=pi_gerstein,scavenge,general

date
echo "Running bash script..."

python run.py

echo "Done."
date

