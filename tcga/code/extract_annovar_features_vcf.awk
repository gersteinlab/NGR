#!/bin/gawk -f

# variant type is given as an argument (default is empty string dealt with as somatic_MC3)
# cancer_tpye is given as an argument (default is empty string dealt with as pancancer, i.e. encompasses all variants in a file)

# returns the value after "=" in an annovar feature"s annotation
# example: feature_annotation = SIFT_score=0.1 -> 0.1
function extract_feature_value(feature_annotation){
	start_index = index(feature_annotation, "=")
	value = substr(feature_annotation, start_index+1, length(feature_annotation))
	return value
}
# updates dictionaries of numeric values
function update_numeric_feature_dict(feature_dict, gene_name, features, feature_index){
	feaure_value = extract_feature_value(features[feature_index])
	if(feaure_value != "." && feaure_value != "")
		feature_dict[gene_name] = feaure_value ";" feature_dict[gene_name]
}

# aggregates a numeric feature"s dictionary results
# example entry and its averaged results: 
# feature_dict["MYC"] = 0.673;0.672;0.671; ->	MYC	0.672	3
function aggregate_results(feature_dict, feature_name){
	feature_str = "gene," feature_name "_avg,N\n" # header
	for(gene in feature_dict){
		split(feature_dict[gene], feature_values, ";")

		gene_total_score = 0
		for(j=1; j<length(feature_values); j++)
			gene_total_score = gene_total_score + feature_values[j]
		
		gene_avg_score = (gene_total_score / (length(feature_values)-1))
		gene_score_num = length(feature_values)-1
		feature_str = feature_str gene "," gene_avg_score ","  gene_score_num "\n"
	}

	return feature_str
}

# prints a categorical feature"s dictionary results
function print_dict_results(feature_dict, feature_name){
	feature_str = "gene," feature_name "_N\n" # header
	for(gene in feature_dict)
		feature_str = feature_str gene ","  feature_dict[gene] "\n"

	return feature_str
}

# checks if the TCGA barcodes correspond to samples of the given cancer type
function check_data_point(data_point, variant_type, cancer_type_barcodes, data_point_functional_reference){
	if((data_point_functional_reference ~ "exonic" || data_point_functional_reference ~ "splicing") && !(data_point_functional_reference ~ "nc")){ # choose coding exonic and splicing variants only		
		if(variant_type == "germline"){ # last columns in germline annotations are formatted differently from somatic ones
			for(i=NF; i>0; i--){
				if($i ~ "TCGA"){
					split($i, feature_vals, ":")
					TCGA_barcode_prefix = substr(feature_vals[length(feature_vals)], 1, 12)
					TCGA_sample_type = substr(feature_vals[length(feature_vals)], 14, 2) # sample types: [0-9]: tumor types, [10-19] normal types, [20+] control types

					if((cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)) && (int(TCGA_sample_type) > 9 && int(TCGA_sample_type) < 15)) # 'normal' samples
						return 1
	
				} else{
					break
				}
			}
		} else{
			split(data_point, feature_vals, "=")
			for(j=12; j<=length(feature_vals); j++){ # TCGA barcodes start at index 12
				if(feature_vals[j] ~ "TCGA-"){
					TCGA_barcode_prefix = substr(feature_vals[j], 1, 12) # unique identifier prefix of a sample in a TCGA barcode
					TCGA_sample_type = substr(feature_vals[j], 14, 2) # sample types: [0-9]: tumor types, [10-19] normal types, [20+] control types
	
					if((cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)) && int(TCGA_sample_type) < 10) # tumor samples
						return 1
				} else{
					break
				}
			}
		}
	} 
		
	return 0
}

BEGIN{
	cancer_type_barcodes = ""
}

{
	if(pass == "barcodes"){ # barcodes file
		# clinical data file (PanCan_ClinicalData_V4_wAIM.txt) with barcodes in the second column
		sample_barcode = $1
		sample_cancer_type = $2
		if(sample_cancer_type == cancer_type)
			cancer_type_barcodes = cancer_type_barcodes sample_barcode "-"
	} else if(NR > 1 && index($1, "#") != 1){ # annotated file
		# split features in annovar annotations
		split($8, features, ";");
			
		annovar_start_index = 2 # in somatic files, annovar_start_index is fixed = 2
		if(variant_type == "germline"){ # in germline, it varies
			for(i=1; i<=length(features); i++){
				if(features[i] ~ "ANNOVAR_DATE"){
					annovar_start_index=i
					break
				}
			}
		}
		
		data_point_functional_reference = substr(features[annovar_start_index+1], 14, length(features[annovar_start_index+1])) # used to select specific functional references (ANNOVAR's Func.refGene) only (e.g. exonic, splicing, etc.)
		include_data_point = check_data_point($10, variant_type, cancer_type_barcodes, data_point_functional_reference) # check if data point corresponds to cancer type under study
		
		if(include_data_point == 1){
			# some genes values are composite; split them
			genes_str = extract_feature_value(features[annovar_start_index+2])
			split(genes_str, gene_names, "\\\\x3b")
			
			for(i=1; i<= length(gene_names); i++){
				# get gene name
				gene_name =  gene_names[i]		
				
				# I. Categorical features
				
				# Exonic function value: synonymous_SNV, nonsynonymous_SNV, stopgain, frameshift_substitution 
				exonic_func_value = extract_feature_value(features[annovar_start_index+4])
				if(exonic_func_value == "nonsynonymous_SNV")
					nonsynonymous[gene_name] = nonsynonymous[gene_name] + 1
				else if(exonic_func_value == "synonymous_SNV")
					synonymous[gene_name] = synonymous[gene_name] + 1
				else if(exonic_func_value == "stopgain")
					stopgain[gene_name] = stopgain[gene_name] + 1
				else if(exonic_func_value == "frameshift_substitution")
					frameshift_substitution[gene_name] = frameshift_substitution[gene_name] + 1
				
				# II. Numeric features
				
				update_numeric_feature_dict(SIFT_scores, gene_name, features, annovar_start_index+8)
				update_numeric_feature_dict(FATHMM_scores, gene_name, features, annovar_start_index+20)
				update_numeric_feature_dict(PROVEAN_scores, gene_name, features, annovar_start_index+23)
				update_numeric_feature_dict(METASVM_scores, gene_name, features, annovar_start_index+26)
				update_numeric_feature_dict(METALR_scores, gene_name, features, annovar_start_index+29)
				update_numeric_feature_dict(MCAP_scores, gene_name, features, annovar_start_index+31+(variant_type == "germline"))	# MCAP scores are 31st annovar features in somatic files, 32nd in germline ones.
				update_numeric_feature_dict(fitCons_scores, gene_name, features, annovar_start_index+45)
				update_numeric_feature_dict(GERP_scores, gene_name, features, annovar_start_index+48)
				update_numeric_feature_dict(PhyloP_scores, gene_name, features, annovar_start_index+52)
				update_numeric_feature_dict(PhyloPhen2_scores, gene_name, features, annovar_start_index+151)
				update_numeric_feature_dict(Eigen_scores, gene_name, features, annovar_start_index+207)
			}
		}
	
		if(FNR % 50000 == 0)
			print "Line " FNR " processed."
	}
}


END{
	# I. Categorical features
	
	print print_dict_results(nonsynonymous, "nonsynonymous") > output_prefix"_nonsynonymous_results.csv"
	print print_dict_results(synonymous, "synonymous") > output_prefix"_synonymous_results.csv"
	print print_dict_results(stopgain, "stopgain") > output_prefix"_stopgain_results.csv"
	print print_dict_results(frameshift_substitution, "frameshift_substitution") > output_prefix"_frameshift_substitution_results.csv"
	
	# II. Numeric features
	
	print aggregate_results(SIFT_scores, "SIFT") > output_prefix"_SIFT_results.csv"
	print aggregate_results(FATHMM_scores, "FATHMM") > output_prefix"_FATHMM_results.csv"
	print aggregate_results(PROVEAN_scores, "PROVEAN") > output_prefix"_PROVEAN_results.csv"
	print aggregate_results(METASVM_scores, "METASVM") > output_prefix"_METASVM_results.csv"
	print aggregate_results(METALR_scores, "METALR") > output_prefix"_METALR_results.csv"
	print aggregate_results(MCAP_scores, "MCAP") > output_prefix"_MCAP_results.csv"
	print aggregate_results(fitCons_scores, "fitCons") > output_prefix"_fitCons_results.csv"
	print aggregate_results(GERP_scores, "GERP") > output_prefix"_GERP_results.csv"
	print aggregate_results(PhyloP_scores, "PhyloP") > output_prefix"_PhyloP_results.csv"
	print aggregate_results(PhyloPhen2_scores, "PhyloPhen2") > output_prefix"_PhyloPhen2_results.csv"
	print aggregate_results(Eigen_scores, "Eigen") > output_prefix"_Eigen_results.csv"
}
