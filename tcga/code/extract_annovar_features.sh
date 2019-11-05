#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=50G
#SBATCH --time=6-23:59:00

date

declare -a cancer_types=('GBM' 'HNSC' 'KICH' 'KIRC' 'KIRP' 'LAML' 'LGG' 'OV' 'PAAD' 'PCPG' 'READ' 'SARC' 'SKCM' 'STAD' 'TGCT' 'THCA' 'THYM' 'UCEC' 'UCS' 'UVM');
declare -a variant_types=('germline');
declare -a txt_files=('germline/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_rmLowQual.vcf.hg19_multianno.txt');
declare -a vcf_files=('germline/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_rmLowQual.vcf.hg19_multianno.vcf');
clinical_file=~/hm444/NGR/tcga/clinical_data/PanCan_ClinicalData_V4_wAIM.txt

for (( j=0; j<${#variant_types[@]}; j++ ));
do
	variant_type=${variant_types[$j]}
	txt_input_file=${txt_files[$j]}
	vcf_input_file=${vcf_files[$j]}

	echo "Variant Type: $variant_type"
	echo "TXT Input File: $txt_input_file"
	echo "VCF Input File: $vcf_input_file"

	for (( i=0; i<${#cancer_types[@]}; i++ ));
	do
		cancer_type=${cancer_types[$i]}
		output_prefix=$variant_type"_"$cancer_type
		echo "Output Prefix: $output_prefix"

		echo "Working on $cancer_type, somatic variants..."
		echo "TXT annotations..."
		awk -v variant_type=$variant_type -v output_prefix=$output_prefix -v cancer_type=$cancer_type -f extract_annovar_features_txt.awk pass="barcodes" $clinical_file pass="" $txt_input_file # to extract number of (1) variants and (2) samples per gene
		date
		
		awk -v variant_type=$variant_type -v output_prefix=$output_prefix -v cancer_type=$cancer_type -f extract_annovar_features_vcf.awk pass="barcodes" $clinical_file pass="" $vcf_input_file # to extract categorical and numeric features
		echo "VCF annotations..."
		date
	done
done

echo "Done."
