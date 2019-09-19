#!/bin/bash

#SBATCH --partition=pi_gerstein
#SBATCH --mem=20G
#SBATCH --time=48:00:00

date

awk '{
        #variant record
        variant = "" $5 "\t" $6 "\t" $7 "\t" $12 "\t" $13 "\t.\t.\t.\t";
        format_col = "GT:";

        # INFO col wil all values separated by "=" as a special character
        info_col = $1 "=" $2 "=" $3 "=" $4 "=" $8 "=" $9 "=" $10 "=" $11;
        for(i=12; i<NF; i++){
                info_col = info_col $i "=";
        }

        info_col = info_col $NF;

        if(NR == 1){ # header
                variant = "#" variant;
                format_col = format_col info_col;
                info_col = "values_of_prev_col_format";
        } else
                format_col = format_col "feature_vals";

        variant = variant format_col "\t" info_col;
        print variant;
}' mc3.v0.2.8.PUBLIC.maf > mc3.v0.2.8.PUBLIC.vcf

date
