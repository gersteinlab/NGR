#!/bin/awk -f

# Filtering PPI networks to only include genes from the scored ones:

{
	if(pass == 1){ # first file: list of genes and scores: e.g. TERT,0.36
		split($0, values, ",")
		gene[values[1]] = values[2]
	} else if (pass == 2){ # second file: PPI as a list of edges: 
		if(gene[$1] != "" && gene[$2] != ""){
			G[$1][$2] = $3 # update the graph
			selected_genes[$1] = 1 # genes to keep from the score list
			selected_genes[$2] = 1
		}
	}
}

END{
	matrix_str = ""
	n = asorti(selected_genes)
	n = 10
	
	# build header
	for(h=1; h<=n; h++){
		matrix_str = matrix_str selected_genes[h]
		
		if(h != n)
			matrix_str = matrix_str "\t"
		else
			matrix_str = matrix_str "\n"
	}
	
	# build matrix
	for(i=1; i<=n; i++){
		gi = selected_genes[i]
		matrix_str = matrix_str gi "\t"
		
		for(j=1; j<=n; j++){
			gj = selected_genes[j]
			
			edge_value = 0
			if(G[gi][gj] != "") # undirected networks
				edge_value = G[gi][gj]
			else if(G[gj][gi] != "")
				edge_value = G[gj][gi]

			matrix_str = matrix_str edge_value
			
			if(j != n)
				matrix_str = matrix_str "\t"
		}
		
		if(i != n)
			matrix_str = matrix_str "\n"
	}
	
	print matrix_str > file_prefix"_filtered_matrix.txt"
}
