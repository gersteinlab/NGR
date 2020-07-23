# To generate gold standard score list files, locally run:
code/generate_gold_standard_lists.Rmd. Check parameters and values inside the Rmd file if needed.

NOTE: current gold standard lists are statistically significant w.r.t. ranking of known genes (last code 
chunk of generate_gold_standard_lists.Rmd). If they are to be used, go fast over the code for a refresher
and sanity check. Only if they do not fit well after experimenting with NGR, you might want to adjust the 
formulas of combining DepMap and CancerMine scores in case 'better' lists are possible and needed. If 
cancer_tpye-specific lists are needed, generate them accordingly and validate with cancer_type-specific 
COSMIC genes per the wilcoxon test in the last code chunk of generate_gold_standard_lists.Rmd. 

# To generate indicator variable feature score files (i.e. indicator variables), locally run the commands below in known_genes/:
awk -F "," 'BEGIN{print "gene,list_membership_COSMIC_v90”} {print$1",1"}' COSMIC_v90_Census_all_Oct_1_19.csv > COSMIC_v90_membership.csv
awk -F "," 'BEGIN{print "gene,list_membership_Vo”} {if(NR>2 && index($1, "classif") == 0 && index($1, "*") == 0){print$1",1"}}' Vogelstein_et_al_list_Jul_15_19_S2A.csv > Vogelstein_list_membership.csv
awk -F "," '{if(NR>2 && index($1, "classif") == 0 && index($1, "*") == 0){print$1",1"}}' Vogelstein_et_al_list_Jul_15_19_S2B.csv >> Vogelstein_list_membership.csv
awk 'BEGIN{print "gene,list_membership_Fr”} {for(i=1; i<=NF; i++){print $i",1"}}' Frampton_et_al_list_Oct_1_19.csv > Frampton_list_membership.csv

# To generate cancer-specific COSMIC lists, 
# locally run commands as per COSMIC_cancer_type_lists_commands.txt