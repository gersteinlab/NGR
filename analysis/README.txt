# To generate final list details used in the final Google sheet, run mobility_lists_analysis.Rmd

# To generate preliminary heatmap and circos data, run calculate_ct_gene_matrix.Rmd

# To generate the heatmaps of enrichment analysis, run generate_gene_lists.Rmd then generate_enrichment_heatmaps.Rmd

# To generate the (final mobility list + driver or COSMIC) gene membership (+ initial score before diffusion) matrix used in creating ppi analysis figure, run generate_gene_lists.Rmd

# To generate other methods lists (Hierarchical HotNet and Katz), run generate_other_methods_lists.pynb

# To compare NGR lists to other methods lists (generated and manually curated ones), run compare_with_other_methods_lists.Rmd


# Cytoscape figure notes/steps:

In Cytoscape:

Enrichment analysis figure:
==


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
- Load HumanNetv2_converted.txt networks (File > Import > Network from File), name columns as Node1, Node2, and edge_score)
- Filter out genes neither in Bailey et al. nor in upward mobility gene list (i. add Bailey_driver_COAD_and_READ_COSMIC_and_upward_mobility_genes_all.txt as attribute from table and rename column to 'include', ii. use 'Select' panel to remove ones without value = 1 in the 'include' column as imported in step i.), ii. cte nodes selected in step ii. 
- Add new columns with the name of the PPI as value; to be used when building composite edges after the merge
- Repeat the three previous steps for STRING_converted.txt
- Scale STRING edge scores to from 0-999 to 0-0.99 (edge_score renamed to raw_edge_score, new column names edge_score); edge_score in HumanNetv2 corresponds to percentile
- Merge networks: Tools > Merge (+ select edge_score and edge_type as lists)
- Import the membership matrix as a table as a table in Nodes view
- Import Merged_Network_Union_node_metrics.csv as a table in Nodes view
- Add a new column for node labels called n_chars =LEN(${name})
- Finalize style of the merged network (edge thickness = average of ppi edge scores, edge color = STRING, HumanNet or Both)
- Duplicate and create a the network according to each cancer type (new columns for node color = driver/COSMIC, upward mobility, or Both)

NOTE: ALL FIGURES ARE GENERATED BASED ON FINAL MOBILITY LISTS, i.e. after filtering w.r.t. 50% threshold in CRISPR and RNAi. Priotized lists are one that include ones without filtering. Whole lists contain mobility lists of all genes in the PPI. Therefore: Whole Lists > Prioritized Lists (in the top chunk after diffusion converges) > Final Lists (after filtering w.r.t. 50% DepMap threshold).

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

PPI analysis figure:
==
Node size: ct_specific initial score = 0.0019-0.3+ mapped to 39.5-431.49, Border width: Degree = 17-773 mapped to 4.0-30.0
Large-small driver node split is based on initial score: < 0.05 = small, 0.05+ = large
More-less connected diver node split is based on degree: < 150 = less connected, 150+ more connected