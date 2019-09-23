#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=50G
#SBATCH --time=6-23:59:00

date

txt_input_file=$1
vcf_input_file=$2
variant_type=$3
output_prefix=$4

echo "TXT Input File: $txt_input_file"
echo "VCF Input File: $vcf_input_file"
echo "Variant Type: $variant_type"
echo "Output Prefix: $output_prefix"

echo "Extracting features from txt file..."
awk -v variant_type=$variant_type -v output_prefix=$output_prefix -f extract_annovar_features_txt.awk $txt_input_file # to extract number of (1) variants and (2) samples per gene
echo "Extraction from txt file done."
date


echo "Extracting features from vcf file..."
awk -v variant_type=$variant_type -v output_prefix=$output_prefix -f extract_annovar_features_vcf.awk $vcf_input_file # to extract categorical and numeric features
echo "Extraction from vcf file done."

echo "Done."
date

