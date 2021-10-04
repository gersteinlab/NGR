#!/bin/bash

#SBATCH --partition=pi_gerstein,general
#SBATCH --job-name=ngr_all_cancer_types
#SBATCH --mem=55GB
#SBATCH --time=6-23:59:59

declare -a cancer_types=('BRCA' 'CESC' 'CHOL' 'COAD' 'ESCA' 'HNSC' 'KICH' 'KIRC' 'KIRP' 'LIHC' 'LUAD' 'LUSC' 'READ' 'PRAD' 'STAD' 'THCA' 'UCEC');
declare -a ppi_networks=('STRING' 'HumanNetv2');

date

for (( j=0; j<${#cancer_types[@]}; j++ ));
do
	cancer_type=${cancer_types[$j]}
	for (( i=0; i<${#ppi_networks[@]}; i++ ));
	do
		ppi_network=${ppi_networks[$i]}
		echo "Running NGR on $cancer_type and $ppi_network..."
		
		python run.py --ppi $ppi_network --ct $cancer_type
	done
done

echo "Done."
date