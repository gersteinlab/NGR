# Select STRING information channels and calculate combined score (selected channels: all except text mining and database in their direct and homology-based forms):
awk 'BEGIN{prior=0.041; feature_inds="1,2,3,4,5,6,7,8,9,10,11"; split(feature_inds, features, ",")} {combined_score = 1.0; interaction_str=""; for(i=1; i<=length(features); i++){if(i>2 && $features[i] > 0) {feature_val = (($features[i]/1000)-prior)/(1-prior); combined_score = combined_score * (1-feature_val)}; interaction_str = interaction_str $features[i] " "} if(FNR == 1){interaction_str = interaction_str "combined_score"} else{combined_score = 1 - combined_score; if(combined_score > 0){combined_score = combined_score + (prior * (1-combined_score)); combined_score = combined_score*1000}; interaction_str = interaction_str combined_score}; if(combined_score > 0){print interaction_str}}' 9606.protein.links.full.v11.0.txt > 9606.protein.links.full.notextmining.nodatabase.v11.0.txt

# Sort STRING edges (optional):
cat 9606.protein.links.full.notextmining.nodatabase.v11.0.txt | awk 'NR<2{print $0;next}{print $0| "sort -k 12,12nr"}' > temp.txt
mv temp.txt 9606.protein.links.full.notextmining.nodatabase.v11.0.txt

# Select high quality edges in STRING (combined_score > 800):
awk '{if($12 > 800 || FNR == 1){print $0}}' 9606.protein.links.full.notextmining.nodatabase.v11.0.txt > 9606.protein.links.full.notextmining.nodatabase.edgescore.800upwards.v11.0.txt

# Convert resulting STRING network into a edgelist file (with 3 columns)
awk '{print $1 " " $2 " " $12}' STRING_9606.protein.links.full.notextmining.nodatabase.edgescore.800upwards.v11.0.txt > STRING_9606.protein.links.notextmining.nodatabase.edgescore.800upwards.v11.0.txt

# Select high quality edges (top 10%) in HumanNetv2:
awk '{if(FNR < 37151){print $0}}' HumanNetv2-FunctionalGeneNetwork_\[FN\].tsv > HumanNetv2-FunctionalGeneNetwork_[FN]_top_10perc_edges.tsv

# add uniform edge weight to HuRI.tsv
awk '{print $1 "\t" $2 "\t1"}' HuRI.tsv > temp.tsv
mv temp.tsv HuRI.tsv

# Replace gene IDs to gene symbols in STRING and HumanNetv2:
Run ppi/code/convert_ppi_IDs.Rmd on macbook twice (with appropriate network index and filenames set)

# Generate HPRD network with weighted edge values (which has gene symbols):
awk -F "\t" 'function ceil(v){ if(v == int(v)){return v} else{return int(v)+1}} {split($8, papers, ",") ;print $1 " " $4 " " ceil(length(papers)/5)/3}' HPRD_BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt > HPRD_converted.txt

# (Optional) Sort edgelists (after moving STRING and HumanNetv2 files from macbook, files with gene symbols end with _converted.txt and are placed in ppi/code/results/)
sort -k 3nr,12 STRING_converted.txt -o STRING_converted.txt
sort -k 3nr,12 HumanNetv2_converted.txt -o HumanNetv2_converted.txt
sort -k 3nr,12 HPRD_converted.txt -o HPRD_converted.txt
sort -k 3nr,12 HuRI_converted.txt -o HuRI_converted.txt

# Remove duplicate edges from STRING (STRING is undirected but originally includes each edge twice as A B score and B A score):
awk '{if($1 < $2){print $0}}' STRING_converted.txt > temp.txt
mv temp.txt > STRING_converted.txt

# Generate matrices from edgelists for the 4 networks (after moving edgelist files from macbook to Farnam):
sbatch convert_ppi_network_to_matrix_script.sh

# Important notes:
Currently, we consider all 4 networks undirected and generate the matrix accordingly.
From literature, it seems STRING and HPRD are undirected networks, HumanNetv2 is directed (as in disease association paper on HumanNet website, but 2013, Nature Method paper by Ideker lab leverages the network as undirected, for example). We consider HuRI undirected.
An edge list includes each edge once: for undirected graphs, an edge is considered twice in the matrix; for directed graphs, once. For current HumanNetv2, there is only few (<10 out of 37,000) edges listed twice in both directions (A B score 1, B A score 2) with same or close scores; we use the max value and deal with each edge as undirected when creating the matrix using the awk script with "undirected" as network type.
