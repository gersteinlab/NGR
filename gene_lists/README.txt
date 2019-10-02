# To generate dependency score files, locally run code/generate_gene_dependency_lists.Rmd locally.

# To generate indicator variable feature score files (i.e. indicator variables), locally run the commands below in known_genes/:
awk -F "," 'BEGIN{print "gene,_COSMIC_v90_membership"} {print$1",1"}' COSMIC_v90_Census_all_Oct_1_19.csv > COSMIC_v90_membership.csv
awk -F "," 'BEGIN{print "gene,Vo_list_membership"} {if(NR>2 && index($1, "classif") == 0 && index($1, "*") == 0){print$1",1"}}' Vogelstein_et_al_list_Jul_15_19_S2A.csv > Vogelstein_list_membership.csv
awk -F "," '{if(NR>2 && index($1, "classif") == 0 && index($1, "*") == 0){print$1",1"}}' Vogelstein_et_al_list_Jul_15_19_S2B.csv >> Vogelstein_list_membership.csv
awk 'BEGIN{print "gene,Fr_list_membership"} {for(i=1; i<=NF; i++){print $i",1"}}' Frampton_et_al_list_Oct_1_19.csv > Frampton_list_membership.csv
