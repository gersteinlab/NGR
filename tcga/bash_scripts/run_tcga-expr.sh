#!/bin/bash

#SBATCH --partition=general
#SBATCH --mem=20GB
#SBATCH --time=29-23:59:00

date

~/Rscript ../code/tcga-diff_exp.R
echo "[Bash script] Rscript Done."

mv *.RData ../code
rm Rplots.pdf

echo "[Bash script] Files moved and deleted as necessary."

echo "[Bash script] Done."
date
