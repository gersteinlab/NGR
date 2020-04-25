# All lines of this README file, commented and un-commented, include important details. Uncommented lines are commands.
# For backup purposes, all scripts referred below are symbolic links. Actuals sript are in ~/hm444/NGR/tcga/code.
# Unless noted otherwise, all sbatch scripts are run locally with start directory being tcga/cortex_data/variant_data/

# To download variant files, run gdc-client download -m manifest_file -t token_file

# MC3 data download:
gdc-client download -m ../manifest_files/MC3_controlled_manifest.txt -t ../gdc-user-token.2019-07-24T19_35_18.053Z.txt
mv MC3 somatic_MC3

# Germline data by Huang et al. download:
gdc-client download -m manifest_files/Germline_Huang_et_al_controlled_manifest.txt -t gdc-user-token.2019-07-24T19_35_18.053Z.txt

# Germline data filtering to remove LowQual variants and aggregate TCGA IDs to reduce file size from 2.7TB to ~46GB
cd germline
mv ../0fbc9ce6-aed8-4bc6-a24f-9cf8229654a1/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz.
rm -r ../0fbc9ce6-aed8-4bc6-a24f-9cf8229654a1
sbatch gunzip_germline_data.sh
sbatch process_germline_variants.sh
rm PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf # remove file (2.7 TB) after processing.

# Convert MC3 public MAF (somatic) to a VCF format accepted by ANNOVAR
cd somatic_MC3/
sbatch convert_mc3_maf_to_vcf.sh

# Download ANNOVAR databases
mkdir annovar
cd annovar
sbatch download_annovar_dbs.sh

# Generate ANNOVAR annotations (run in tcga/cortex_data/variant_data/)
sbatch generate_annovar_annotations.sh germline/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_rmLowQual.vcf
sbatch generate_annovar_annotations.sh somatic_MC3/mc3.v0.2.8.PUBLIC.vcf

## Aggregate features
# Extract features on Pancancer and tissue-specific levels from ANNOVAR annotations (run in tcga/cortex_data/variant_data/)
# Note: this command takes 14+ days to extract feature for each cancer type. For faster processing, modify extract_annovar_features.sh to run in paralell.
sbatch extract_annovar_features.sh
rename __ _ *
mv *results.csv annovar/results/aggregate_results/

# Normalize variant-related features by gene length:
# GRCh37_bioMart_bioConductor_package_all_121819.txt reference file used inside the shell script was generated using tcga/code/extract_gene_lengths.Rmd
sbatch normalize_varscores_by_gene_length.sh

## Sample-level Features
# Notes: We focus when generating patient-level matrices here on coding variants, exonic or splicing only.

# generate genomics  matrices (sample-level, gene)
sbatch extract_annovar_matrices_txt.sh
rename __ _ *
mv *results.csv annovar/results/matrix_results

# transform matrices (normalization by gene length + discretization per 1-4 scale)
Interactively on mac, run: tcga/code/normalize_matrix_by_gene_length.Rmd

