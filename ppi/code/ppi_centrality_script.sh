#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --job-name=ppi_centrality
#SBATCH --cpus-per-task=4
#SBATCH --time=13-23:59:00
#SBATCH --mem=50GB

start_time=$(date +%s)
date

echo "$@"
date
python ppi_centrality.py "$@"
date

end_time=$(date +%s)

hours=$(( ($end_time-$start_time) / ( 3600 )))
minutes=$(( ($end_time-$start_time) / ( 60 )))
seconds=$(( ($end_time-$start_time) % 60 ))

echo "Execution time: " $hours "h" $(($minutes % (60))) "m" $seconds "s"
