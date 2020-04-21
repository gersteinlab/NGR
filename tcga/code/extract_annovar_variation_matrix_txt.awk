#!/bin/gawk -f

BEGIN{
	cancer_type_barcodes = ""
}

{
	if(pass == "barcodes"){ # barcodes file: build a string of cancer type-specific barcode
		# clinical data file (PanCan_ClinicalData_V4_*.txt) with barcodes in the second column
		sample_barcode = $1
		sample_cancer_type = $2

		if(sample_cancer_type == cancer_type)
			cancer_type_barcodes = cancer_type_barcodes sample_barcode "-"
	} else if(FNR > 1){ # annotated annovar txt file beyond header
		if(($6 ~ "exonic" || $6 ~ "splicing") && !($6 ~ "nc")){
			split($7, gene_names, ";") # to account for certain records where multiple genes are separated by semi-colons
			for(i in gene_names){
				gene = gene_names[i]
				
				# samplepergene
				if(variant_type == "germline"){ # last columns in germline annotations are formatted differently from somatic ones
					for(i=NF; i>0; i--){
						if($i ~ "TCGA"){
							split($i, feature_vals, ":")
	
							TCGA_barcode = substr(feature_vals[length(feature_vals)], 1, 16) # barcode is a sample ID
							TCGA_barcode_prefix = substr(TCGA_barcode, 1, 12) # barcode_prefix is a patient ID (a patient might have multiple samples)
							TCGA_sample_type = substr(TCGA_barcode, 14, 2) # sample types: [0-9]: tumor types, [10-19] normal types, [20+] control types
							
							if((cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)) && (int(TCGA_sample_type) > 9 && int(TCGA_sample_type) < 15)){ # 'normal' samples
								gene_list[gene] = 1 # for output generation purposes
								TCGA_barcode_list[TCGA_barcode] = 1 # for output generation purposes
								
								if(variant_count[TCGA_barcode][gene] == "")
									variant_count[TCGA_barcode][gene] = 1
								else
									variant_count[TCGA_barcode][gene] = variant_count[TCGA_barcode][gene] + 1
							}
						} else{
							break
						}
					}
				} else{
					split($NF, feature_vals, "=")
					for(j=12; j<=length(feature_vals); j++){ # TCGA IDs start at index 12
						if(feature_vals[j] ~ "TCGA-"){												
							TCGA_barcode = substr(feature_vals[j], 1, 16) # barcode is a sample ID
							TCGA_barcode_prefix = substr(TCGA_barcode, 1, 12) # barcode_prefix is a patient ID (a patient might have multiple samples)
							TCGA_sample_type = substr(TCGA_barcode, 14, 2) # sample types: [0-9]: tumor types, [10-19] normal types, [20+] control types
							
							if((cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)) && int(TCGA_sample_type) < 10){ # tumor samples
								gene_list[gene] = 1 # for output generation purposes
								TCGA_barcode_list[TCGA_barcode] = 1 # for output generation purposes
								
								if(variant_count[TCGA_barcode][gene] == "")
									variant_count[TCGA_barcode][gene] = 1
								else
									variant_count[TCGA_barcode][gene] = variant_count[TCGA_barcode][gene] + 1
							}
						} else{
							break
						}
					}
				}
			}
	
			if(FNR % 50000 == 0)
				print "Line " FNR " processed."
		}
	}
}

END{
	# generate three files: barcode index (order of patient IDs in rows), gene index (order of columns), matrix (matrix values, numerics only)

	# barcode index
	n_b = asorti(TCGA_barcode_list) # sort by values of indices (i.e. barcodes)
	barcode_index_str = ""
	for(b=1; b<=n_b; b++)
		barcode_index_str = barcode_index_str TCGA_barcode_list[b] "\n"
		
	print barcode_index_str > output_prefix"_matrix_barcode_index.csv"
	print "Barcode index of " output_prefix " created."
		
	# gene index
	n_g = asorti(gene_list)
	gene_index_str = ""
	for(g=1; g<= n_g; g++)
		gene_index_str = gene_index_str gene_list[g] "\n"
	
	print gene_index_str > output_prefix"_matrix_gene_index.csv"
	print "Gene index of " output_prefix " created."

	# matrix
	matrix_str = ""
	for(b=1; b<=n_b; b++){
		barcode = TCGA_barcode_list[b]

		for(g=1; g<= n_g; g++){
			gene = gene_list[g]
			
			if(variant_count[barcode][gene] > 0)
				matrix_str = matrix_str variant_count[barcode][gene] " "
			else
				matrix_str = matrix_str "0 "
		}
		
		matrix_str = matrix_str "\n"
	}	
	
	print matrix_str > output_prefix"_matrix.txt"
	print "Matrix of " output_prefix " created."
}
