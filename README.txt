All code is on GitHub.

All data sufficient to run the code is on Farnam. In other words, Farnam is sufficient to reproduce all results.

More data/analyses might be found locally (e.g. gene_gene_matrix for more pathway databases than selected one (KEGG))

Feature extraction is done within each subfolder. Each subfolder has its code/ directory.

Each directory corresponds to a directory/data type. Scripts are to be run locally unless noted otherwise.

pathways/:
To extract and analyze pathways: run code/extract_pathways.Rmd

tcga/:
To fetch clinical data: run code/tcga-clinical.Rmd. Current files in clinical_data/ are manually created symbolic links.
To generate differential expression information: run bash_scripts/run_tcga-expr.sh.
Note on differential expression analysis: currently, emphasis is put on FDR<=0.05 only; no logFC threshold. If logFC threshold is to be enforced, update and rerun code/tcga-diff_exp.R to use glmTreat() for model training: see section 2.12 in edgeR manual for more details.

ppi/:
To generate the PPI file with IDs converted to gene symbols, run code/convert_ppi_IDs.Rmd