#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=50Gb
#SBATCH --time=5-23:00:00

date

echo "VarPerGene files..."
for f in annovar/results/aggregate_results/*_varpergene* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

echo "SamplePerGene files..."
for f in annovar/results/aggregate_results/*_samplepergene* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

echo "Synonymous files..."
for f in annovar/results/aggregate_results/*_synonymous* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

echo "Non-Synonymous files..."
for f in annovar/results/aggregate_results/*_nonsynonymous* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

echo "Stop-Gain files files..."
for f in annovar/results/aggregate_results/*_stopgain* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

echo "Frameshift-Sustitution files..."
for f in annovar/results/aggregate_results/*_frameshift_substitution* 
do 
	awk -f normalize_subscores_by_gene_length.awk pass=1 annovar/GRCh37_bioMart_bioConductor_package_all_121819.txt pass=2 $f
	echo $f " normalized successfully."
done

date

echo "Done."
