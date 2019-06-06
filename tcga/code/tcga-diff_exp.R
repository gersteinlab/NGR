display_diff_exp_stats <- function(lrt, p_val_threshold=0.05, 
                                   FDR_threshold=0.05, logFC_threshold=1){
  # summary stats/counts
  summary(decideTests(lrt))
  
  gene_table <- topTags(lrt, n=dim(lrt$table)[1])
  gene_table <- gene_table[order(gene_table$table$FDR), ]
  
  gene_table <- gene_table[abs(gene_table$table$logFC) >= logFC_threshold,] 
  print(paste("# of diff. exp. logFC genes:", dim(gene_table)[1]))

  gene_table <- gene_table[gene_table$table$PValue <= p_val_threshold &
                             abs(gene_table$table$logFC) >= logFC_threshold,]
  
  print(paste("# of diff. exp. logFC genes, filtered by PValue:", dim(gene_table)[1]))
  
  print(paste("# of diff. exp. logFC genes, filtered by FDR (q_value):", dim(gene_table)[1]))
  gene_table <- gene_table[gene_table$table$FDR <= FDR_threshold &
                             abs(gene_table$table$logFC) >= logFC_threshold,] 

  #print(gene_table) # differentially expressed genes with FDR < 0.05 and logFC > 1
  return(gene_table)
}

library(edgeR)

# Tutorial: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf
expression_data_dir <- '../expression_data/'
cancer_types <- c('brca', 'luad', 'lusc', 'prad', 'lihc', 'chol', 'meso', 'dlbc')
#cancer_types <- c('brca')
tumor_exp_files <- c(paste(expression_data_dir, cancer_types, '-rsem-count-tcga-t.txt', sep=''))
healthy_exp_files <- c(paste(expression_data_dir, cancer_types, '-rsem-count-tcga.txt', sep=''))

for(i in 1:length(cancer_types)){
  tryCatch({
    print(paste("Execution for ", cancer_types[i], " started. ", "[", Sys.time(), "]", sep=""))
    
    t_counts <- read.table(tumor_exp_files[i], header=T)
    h_counts <- read.table(healthy_exp_files[i], header=T)
    
    # Hugo (gene) symbol is in col 1, Entrez ID is in col 2, counts are in cols [3, ]
    # Genes must be of the same order in both files. See ../expressiom_data/ for examples
    gene_symbols <- t_counts[, 1]
    gene_entrez_ids <- t_counts[, 2]
    
    t_counts <-  t_counts[, 3:ncol(t_counts)]
    rownames(t_counts) <- gene_symbols
    
    h_counts <-  h_counts[, 3:ncol(h_counts)]
    rownames(h_counts) <- gene_symbols
    
    # Code block that uses only tumor and matched samples
    # TCGA sample ID is in characters [9,12]. Keep samples from participants from whom both healthy and tumor tissues were taken.
    # See TCGA sample code: https://docs.gdc.cancer.gov/Encyclopedia/pages/images/TCGA-TCGAbarcode-080518-1750-4378.pdf
    # intersect_samples <- intersect(substr(colnames(t_counts), 9, 12), substr(colnames(h_counts), 9, 12))
    # h_select_samples <- which(substr(colnames(h_counts), 9, 12) %in% intersect_samples)
    # t_select_samples <- which(substr(colnames(t_counts), 9, 12) %in% intersect_samples)
    # 
    # # double check: keeps one occurrence for healthy and tumor tissue per participant
    # h_select_samples <- match(intersect_samples, substr(colnames(h_counts), 9, 12))
    # t_select_samples <- match(intersect_samples, substr(colnames(t_counts), 9, 12))
    # 
    # h_counts <-  h_counts[, h_select_samples]
    # t_counts <-  t_counts[, t_select_samples]
    # 
    # # order oclumns of healthy and tumor (resulting data frame has both samples and genes in the same order)
    # h_counts <- h_counts[, order(substr(colnames(h_counts), 9, 12))]
    # t_counts <- t_counts[, order(substr(colnames(t_counts), 9, 12))]
    #
    # merge healthy and tumor, create DGEList object
    #counts <- cbind(h_counts, t_counts)[, order(c(seq(ncol(h_counts)), seq(ncol(t_counts))))] # alternate columns in cbind (to match designMatrix below)

    counts <- cbind(h_counts, t_counts)
    dgList <- DGEList(counts=counts, genes=gene_symbols)
    
    print(paste("dgList object created. ", "[", Sys.time(), "]", sep=""))
    
    # Filtering: retain only those genes that are represented at least 1cpm reads in at least two samples (cpm=counts per million).
    countsPerMillion <- cpm(dgList)
    #summary(countsPerMillion)
    #'summary' is a useful function for exploring numeric data; eg. summary(1:100)
    countCheck <- countsPerMillion > 1
    keep <- which(rowSums(countCheck) >= 2)
    dgList <- dgList[keep, ]
    #summary(cpm(dgList)) #compare this to the original summary
    
    print(paste("Filtering done. ", "[", Sys.time(), "]", sep=""))
    
    dgList <- calcNormFactors(dgList, method="TMM") # normalization step
    save.image(paste(cancer_types[i], "checkpoint1.RData", sep="_"))
    
    sampleType <- factor(c(rep(1, ncol(h_counts)), rep(2, ncol(t_counts)))) # 1 for healthy, 2 for tumor
    designMat <- model.matrix(~sampleType)
    
    dgList <- estimateGLMCommonDisp(dgList, design=designMat) # estimating dispersions
    print(paste("Common dispersions calculated. ", "[", Sys.time(), "]", sep=""))
    save.image(paste(cancer_types[i], "checkpoint2.RData", sep="_"))

    dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
    print(paste("Tagwise dispersions calculated. ", "[", Sys.time(), "]", sep=""))
    save.image(paste(cancer_types[i], "checkpoint3.RData", sep="_"))
    
    fit <- glmQLFit(dgList, designMat)
    lrt <- glmQLFTest(fit, coef=2)
    
    #edgeR_result <- topTags(lrt)
    #save(topTags(lrt,n=15000)$table, file='edgeR_Result.RData')
    
    print(paste("Model fit. ", "[", Sys.time(), "]", sep=""))
    print(lrt$table[order(lrt$table[, "PValue"]), ])
    save.image(paste(cancer_types[i], "checkpoint4.RData", sep="_"))
    
    deGenes <- decideTestsDGE(lrt, p=0.001)
    deGenes <- rownames(lrt)[as.logical(deGenes)]
    plotSmear(lrt, de.tags=deGenes)
    abline(h=c(-1, 1), col=2)
    print(paste("Plot generated. ", "[", Sys.time(), "]", sep=""))
    save.image(paste(cancer_types[i], "checkpoint5.RData", sep="_"))
    
    display_diff_exp_stats(lrt)
  }, error = function(e) {
    print(paste("Error with ", cancer_types[i], sep=""))
    print(e)
  })
}
