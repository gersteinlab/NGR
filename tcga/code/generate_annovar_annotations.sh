#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=50G
#SBATCH --time=6-23:59:00

date

input_file=$1
echo "Input File: $1"

module load ANNOVAR

table_annovar.pl $input_file annovar/ -buildver hg19 -protocol refGene,GDI_full_10282015,genomicSuperDups,dbnsfp35c,snp129,gnomad211_exome,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,1000g2015aug_eas,esp6500siv2_all,popfreq_all_20150413,gme,avsnp150,exac03,clinvar_20190305,cosmic70,regsnpintron,mitimpact24,eigen,revel,dbscsnv11,gwava -operation g,f,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -vcfinput

date
