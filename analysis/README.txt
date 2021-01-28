# To generate final list details used in the final Google sheet, run mobility_lists_analysis.Rmd

# To generate preliminary heatmap and circos data, run calculate_ct_gene_matrix.Rmd

# To generate the heatmaps of enrichment analysis, run generate_gene_lists.Rmd then generate_enrichment_heatmaps.Rmd

# To generate the (final mobility list + driver or COSMIC) gene membership (+ initial score before diffusion) matrix used in creating ppi analysis figure, run generate_gene_lists.Rmd

# To generate other methods lists (Hierarchical HotNet and Katz), run HHotNet/NGR_related/NGR_parallel.sh and generate_other_methods_lists.pynb

# To compare NGR lists to other methods lists (generated and manually curated ones), run compare_with_other_methods_lists.Rmd

# To generate results on popular genes in most connected entities in EnrichmentMap, of the pancancer enrichment analysis, run extract_cluster_popular_UMGs.awk in analysis/Fig_enrichment_analysis_results/pancancer_results/ on macbook  

# Cytoscape figure notes/steps:

In Cytoscape:

Cancer type-specific panel:
=

To generate node metrics (degree, neighborhood connectedness, betweenness centrality, etc.) in merged network (Union of STRING and HumanNetv2):
[Run once, save results in Merged Network_Union_node_metrics.csv]
- Load HumanNetv2_converted.txt network
- Load STRING_converted.txt network
- Merge networks (Union)
- Tools > Network Analyzer > Network Analysis > Analyze network
- Export results (Nodes table) into Merged_Network_Union_node_metrics.csv
- Manually remove SUID, is.selected, matching.attribute, shared name columns and place name column as first column (to be used as key by default)

To generate cancer type specific networks:
- Load HumanNetv2_converted*.txt networks (File > Import > Network from File), name columns as Node1, Node2, and edge_score)
- Filter out genes neither in Bailey et al. nor in upward mobility gene list (i. add Bailey_driver_COAD_and_READ_COSMIC_and_upward_mobility_genes_all.txt as attribute from table and rename column to 'include', ii. use 'Select' panel to remove ones without value = 1 in the 'include' column as imported in step i.), ii. cte nodes selected in step ii. 
- Add new columns with the name of the PPI as value; to be used when building composite edges after the merge
- Repeat the three previous steps for STRING_converted.txt
- Scale STRING edge scores to from 0-999 to 0-0.99 (edge_score renamed to raw_edge_score, new column names edge_score); edge_score in HumanNetv2 corresponds to percentile
- Merge networks: Tools > Merge (+ select edge_score and edge_type as lists)
- Import the membership matrix as a table in Nodes view
- Import Merged_Network_Union_node_metrics.csv as a table in Nodes view
- Add a new column for node labels called n_chars =LEN(${name})
- Finalize style of the merged network (edge thickness = average of ppi edge scores, edge color = STRING, HumanNet or Both)

For each cancer type:
	- Duplicate and create a the network according to each cancer type (by deleting nodes that are neither driver/COSMIC nor upward mobility)
	- Select cancer type-specific UMGs and create a circular layout for UMGs in the middle
	- Remove inter-group edges (within drivers and UMGs) (https://groups.google.com/g/cytoscape-helpdesk/c/BVyq3-ctXrw)
	- Per sepcs below: 
		- Change node color, node size, and label size to correspond to the columns of the cancer type under study
		- Create a new edge_category column with 4 types (Small_thinBorder, Small_thickBorder, Large_thinBorder, Large_thickBorder)
		- Create circular layout for each of the four categories based on splits wrt initial score and degree
		- Place subgroups clockwise and unconnected subgroup top right
		- Double check node width is set to degree as per specs below
		- Check for genes with near 0 score, be it ones with connection in the figure or else

Specs:

Font size: 0-0.3+ cancer type initial score mapped to 20-140
Node size: ct_specific initial score = 0.0019, 0.025, 0.075, 0.1, 0.7+ mapped to 114, 169, 400, 758, 908
Border width: Degree = 0-17 to 773+ (e.g. BRCA: 773+) mapped to 4.0 to 30.0
Large-small driver node split is based on initial score: < 0.075 = small, 0.075+ = large
More-less connected diver node split is based on degree (i.e. border thickness): < 150 = less connected, 150+ more connected
Node colors: #C14215 for drivers of the cancer type, #F5F5DC for upward mobility genes for that type
Edge colors: Olive: #AE9C45, Pink-Orange: #F5793A, Purple: #A78BFB, Dark blue: #000033

NOTE: ALL FIGURES ARE GENERATED BASED ON FINAL MOBILITY LISTS, i.e. after filtering w.r.t. 50% threshold in CRISPR and RNAi. Priotized lists are one that include ones without filtering. Whole lists contain mobility lists of all genes in the PPI. Therefore: Whole Lists > Prioritized Lists (in the top chunk after diffusion converges) > Final Lists (after filtering w.r.t. 50% DepMap threshold).

Enrichment analysis figure:
==

Pancancer panel:
=

Node size: gs_size = 3-20+ mapped to 90-180, Edge width: Overlap size = 2-25+ mapped to 10-45  thickness, Node color: Pathway = #000066, BP = #0F6BDF, MF = #FE7B4E

Custom clusters:
knp: Known Cancer Pathways and Processes
prf: Proliferation
adh: Cell Adhestion and Cellular Matrix
bld: Blood (and Angiogenesis)
trs: Transcription- and Translation-related
bnd: Binding
dvm: Development and Maintenance
imn: Immune System-Related
nrv: Nervous System-Related
cnr: Cancers
dsi: Non-Cancer Diseases and Infections

# Method Comparisons Figure:
Run the second section of compare_with_other_methods_lists.Rmd