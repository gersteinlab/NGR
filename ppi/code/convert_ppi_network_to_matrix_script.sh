#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --job-name=ppi_network_conversion_to_matrix
#SBATCH --cpus-per-task=4
#SBATCH --time=13-23:59:00
#SBATCH --mem=50GB

start_time=$(date +%s)
date

echo "Working on STRING network..."
awk -v file_prefix="STRING_converted" -f convert_ppi_network_to_matrix.awk results/STRING_converted.txt
echo "Working on HumanNet_v2 network...."
awk -v file_prefix="HumanNetv2_converted" -f convert_ppi_network_to_matrix.awk results/HumanNetv2_converted.txt
date

end_time=$(date +%s)

hours=$(( ($end_time-$start_time) / ( 3600 )))
minutes=$(( ($end_time-$start_time) / ( 60 )))
seconds=$(( ($end_time-$start_time) % 60 ))

echo "Execution time: " $hours "h" $(($minutes % (60))) "m" $seconds "s"
