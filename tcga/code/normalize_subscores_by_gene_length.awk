#!/bin/gawk -f

BEGIN{
	output_str = ""
	output_filename = ""
	avg_gene_length = 28679 # calculated manually from the biomart database (see file in pass=1)
}

{
	if(pass == 1){ # biomart file generated using extract_gene_lengths.Rmd (gene symbol in col 2, gene length in col 7)
		gene_length[$2] = $7
	} else if (pass == 2){ # score file to be normalized (comma separated, gene_symbol,N)
		if(FNR == 1){ # header and file name setting
			output_filename = FILENAME
			output_str = output_str "gene,normalized_by_length_orbyavglengthifnotfound_N, normalized_by_length_or0ifnotfound_N,unnormalized_N\n"
		} else if($0 ~ ","){
			split($0, values, ",")
		
			gene_symbol = values[1]
			N = values[2]
	
			if(gene_length[gene_symbol] != "" && gene_symbol != "gene" && length(gene_symbol) >= 2){
				normalized_N = (N/gene_length[gene_symbol])
				normalized_N_ifnotfound = normalized_N
			} else{ # not found in bioMart; divide by average length for first normalized score, replace by 0 for second normalized score
				normalized_N = (N/avg_gene_length)
				normalized_N_ifnotfound = 0
			}
			
			output_str = output_str gene_symbol "," normalized_N "," normalized_N_ifnotfound "," N "\n"
		}
	}
}

END{
	print output_str > output_filename
}

