#!/bin/gawk -f

{
	if(NR > 1){
		# varpergene and samplepergene
		split($7, gene_names, ";") # to account for certain records where multiple genes are separated by semi-colons
		for(i in gene_names){
			gene = gene_names[i]
			
			# varpergene
			if(gene_count[gene] == "")
				gene_count[gene] = 1
			else
				gene_count[gene] = gene_count[gene] + 1
				
			# samplepergene
			if(variant_type == "germline"){ # last columns in germline annotations are formatted differently from somatic ones
				for(i=NF; i>0; i--){
					if($i ~ "TCGA"){
						split($i, feature_vals, ":")
						TCGA_ID = feature_vals[length(feature_vals)]

						if(sample_count[gene][TCGA_ID] == "")
							sample_count[gene][TCGA_ID] = 1
						else
							sample_count[gene][TCGA_ID] = sample_count[gene][TCGA_ID] + 1
					} else{
						break
					}
				}
			} else{
				split($NF, feature_vals, "=")
				for(j=12; j<=NF; j++){ # TCGA IDs start at index 12
					if(feature_vals[j] ~ "TCGA-"){
						TCGA_ID = substr(feature_vals[j], 1, 16) # unique identifier of a sample in a TCGA ID
						
						if(sample_count[gene][TCGA_ID] == "")
							sample_count[gene][TCGA_ID] = 1
						else
							sample_count[gene][TCGA_ID] = sample_count[gene][TCGA_ID] + 1
					} else{
						break
					}
				}
			}
		}
		
		if(FNR % 50000 == 0)
			print "Line " NR " processed."
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
