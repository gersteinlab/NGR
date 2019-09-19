#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=10G
#SBATCH --time=96:00:00

date

awk '{if(index($1, "#") == 1){print $0}; if($1 ~ "#CHROM"){for(i=10; i<=NF; i++){TCGA_ids[i] = $i;}} else if($7 == "PASS" || $7 == "."){TCGA_info=""; 
for(i=1; i<=NF; i++){ 
if(i<=8) {printf $i "\t"}
else if(i == 9){printf $i":sample_name"} # FORMAT columns, all defaults and sample_name
else{if(index($i, ".:.") == 0){printf "\t" $i":"TCGA_ids[i]}}
} 
print "" }}' PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf > PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_rmLowQual.vcf

date
