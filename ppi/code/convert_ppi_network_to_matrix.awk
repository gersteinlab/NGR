#!/bin/awk -f

# Filtering PPI networks to only include genes from the scored ones:

{
	 	# PPI as a list of edges:
		G[$1][$2] = $3 # update the graph
		genes[$1] = 1 # genes to keep from the score list
		genes[$2] = 1
}

END{
	n = asorti(genes)

	# index file (ordered list of genes)
	index_str = ""
	for(h=1; h<=n; h++)
		index_str = index_str genes[h] "\n"
	
	printf index_str > file_prefix"_matrix_index.txt"

	# build matrix
	matrix_str = ""
	for(i=1; i<=n; i++){
		gi = genes[i]		
		for(j=1; j<=n; j++){
			gj = genes[j]
			
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
	
	print matrix_str > file_prefix"_matrix.txt"
}
