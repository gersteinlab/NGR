#!/bin/gawk -f

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
	} else if(FNR > 1){ # annotated file
		#print FNR		
		# varpergene and samplepergene
		split($7, gene_names, ";") # to account for certain records where multiple genes are separated by semi-colons
		for(i in gene_names){
			gene = gene_names[i]

			include_in_gene_count = 0 # flag on whether to include the variant in varpergene feature
			
			# samplepergene
			if(variant_type == "germline"){ # last columns in germline annotations are formatted differently from somatic ones
				for(i=NF; i>0; i--){
					if($i ~ "TCGA"){
						split($i, feature_vals, ":")

						TCGA_barcode = substr(feature_vals[length(feature_vals)], 1, 16)
						TCGA_barcode_prefix = substr(TCGA_barcode, 1, 12)
						
						if(cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)){
							if(sample_count[gene][TCGA_barcode] == "")
								sample_count[gene][TCGA_barcode] = 1
							else
								sample_count[gene][TCGA_barcode] = sample_count[gene][TCGA_barcode] + 1
							
							include_in_gene_count = 1 # update flag for gene count
						}
					} else{
						break
					}
				}
			} else{
				split($NF, feature_vals, "=")
				for(j=12; j<=length(feature_vals); j++){ # TCGA IDs start at index 12
					if(feature_vals[j] ~ "TCGA-"){
						TCGA_barcode = substr(feature_vals[j], 1, 16) # unique identifier of a sample in a TCGA ID
						TCGA_barcode_prefix = substr(TCGA_barcode, 1, 12) 

						if(cancer_type == "" || (cancer_type != "" && cancer_type_barcodes ~ TCGA_barcode_prefix)){
							if(sample_count[gene][TCGA_barcode] == "")
								sample_count[gene][TCGA_barcode] = 1
							else
								sample_count[gene][TCGA_barcode] = sample_count[gene][TCGA_barcode] + 1

							include_in_gene_count = 1 # update flag for gene count
						}
					} else{
						break
					}
				}
			}
			
			# varpergene
			if(include_in_gene_count){
				if(gene_count[gene] == "")
					gene_count[gene] = 1
				else
					gene_count[gene] = gene_count[gene] + 1
			}
		}

		if(FNR % 50000 == 0)
			print "Line " FNR " processed."
	}
}

END{
	# varpergene
	varpergene_feature = "gene,variant_count\n" # number of variants per gene
	for(gene in gene_count)
		varpergene_feature = varpergene_feature gene "," gene_count[gene] "\n"

	print varpergene_feature > output_prefix"_varpergene_results.csv"
	
	#samplepergene
	samplepergene_feature = "gene,sample_count\n" # number of samples per gene
	for(gene in sample_count)
		samplepergene_feature = samplepergene_feature gene "," length(sample_count[gene]) "\n"
		
	print samplepergene_feature > output_prefix"_samplepergene_results.csv"
}
