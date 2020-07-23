All code is on GitHub.

All data sufficient to run the code is on Farnam. In other words, Farnam is sufficient to reproduce all results.

More data/analyses might be found locally (e.g. gene_gene_matrix for more pathway databases than selected one (KEGG))

Feature extraction is done within each subfolder. Each subfolder has its code/ directory.

Each directory corresponds to a directory/data type. Scripts are to be run locally unless noted otherwise.

Pathway extraction and analysis:
To extract and analyze pathways: run extract_pathways.Rmd in pathways/code/

TCGA data and feature extraction:
In tcga/:
Clinical data in clinical_data/ are provided by Tao Qing of Pusztai lab; see clinical_data/REAMDE.txt if needed.
To generate differential expression information: run code/run_tcga-expr.sh locally. Move results and slurm log to code/results/ directory.
Note on differential expression analysis: currently, emphasis is put on FDR<=0.05 only; no logFC threshold. If logFC threshold is to be enforced, update and rerun code/tcga-diff_exp.R to use glmTreat() for model training: see section 2.12 in edgeR manual for more details.
Differential expression analysis results are in code/results/

For data download, annotation (using ANNOVAR), and feature extraction in somatic and germline TCGA variants, see tcga/README.txt

For PPI ID conversion, processing, and convesion to matrices, see ppi/README.txt

(Optional, ~ OBSOLETE as new PPIs have been generated and sample-level results are used) PPI metric (e.g. centrality) generation:
To generate betweenness centrality results, run (full commands in NGR Board sheet) ppi_centrality_script.sh in ppi/code
Betweenness centrality results are in ppi/code/results

For gold standard list generation, see gene_lists/README.txt.

Combined score generation:
To merge features and generate combined scores to be used as inputs to the method, locally run in base/:
Note: This R script is usually run on macbook but should execute successfuly on Farnam as needed if all R packages are installed.
module load R
Rscript merge_features.R -v somatic_MC3
Rscript merge_features.R -v germline

PPI network matrices:
To generate PPI matrices to be used by the Python script of the method, run sbatch convert_ppi_network_to_matrix_script.sh in ppi/code/

For result generation, see method/README.txt

