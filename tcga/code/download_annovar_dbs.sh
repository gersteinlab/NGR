#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=10G
#SBATCH --time=96:00:00

date

module load ANNOVAR

date

annotate_variation.pl -webfrom annovar -downdb refGene -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb knownGene -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb exac03 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb dbnsfp35c -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb gnomad211_exome -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb gnomad211_genome -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb avsnp150 -buildver hg19 .

date

annotate_variation.pl -webfrom annovar -downdb clinvar_20190305 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb cosmic70 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb regsnpintron -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb eigen -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb gwava -buildver hg19 .

date

annotate_variation.pl -webfrom annovar -downdb 1000g2015aug -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb gme -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb revel -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb dbscsnv11 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb esp6500siv2_all -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb snp129 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb popfreq_all_20150413 -buildver hg19 .
annotate_variation.pl -webfrom annovar -downdb mitimpact24 -buildver hg19 .
date

annotate_variation.pl -downdb -buildver hg19 genomicSuperDups .
wget http://www.openbioinformatics.org/annovar/download/GDI_full_10282015.txt.gz
gunzip GDI_full_10282015.txt.gz
mv GDI_full_10282015.txt hg19_GDI_full_10282015.txt

wget http://www.openbioinformatics.org/annovar/download/RVIS_ExAC_4KW.txt.gz
gunzip RVIS_ExAC_4KW.txt.gz
mv RVIS_ExAC_4KW.txt hg19_RVIS_ExAC_4KW.txt

wget http://www.openbioinformatics.org/annovar/download/LoFtool_scores.txt.gz
gunzip LoFtool_scores.txt.gz
mv LoFtool_scores.txt hg19_LoFtool_scores.txt

date

echo "ANNOVAR databases downloaded."
date
