# functions shareable among multiple scripts, some of which already used in mobility_lists_analysis.Rmd and compare_depmap_scores_across_methods.Rmd

select_genes <- function(df, mobility_score_threshold=1500, final_rank_threshold=1000){
  sdf <- df[df$mobility_score >= mobility_score_threshold & df$final_rank <= final_rank_threshold, ]
  sdf <- sdf[order(sdf$final_rank), ]
  return(sdf)
}

# returns known genes in gene_list
get_known_genes <- function(gene_list, known_genes){
  intersecting_known_genes <- intersect(gene_list, known_genes)
  if(length(intersecting_known_genes) > 0)
    return(capture.output(cat(intersecting_known_genes, sep=', ')))
  else
    return('None')
}

# returns a vector with integer values 0, -1 and +1; 
# 0 if gene is not differentially expressed; 1 over-expressed, -1 underexpressed
get_expr_vector <- function(input_genes, expr_dir, cancer_type){
  expr_data <- read.table(paste(expr_dir, paste(tolower(cancer_type), 'expression', 'results.csv', sep='_'), sep=''), header=T, stringsAsFactors=F, sep=',') # DE genes
  de_inds <- match(input_genes, expr_data$gene)
  
  expr_vector <- expr_data[de_inds, paste('logFC', tolower(cancer_type), sep='_')]
  expr_vector[is.na(expr_vector)] <- 0
  expr_vector <- sign(expr_vector)
  
  return(expr_vector)
}

# returns a df with 2 cols for ALL genes: minimum score across all cell lines in RNAi and CRISPR screens
calculate_depmap_min_df_all_genes <- function(depmap_dir){
  depmap_filename <- paste(depmap_dir, 'Achilles_gene_effect_all.csv', sep='')
  crispr_depmap_df_all_genes <- read.table(depmap_filename, header=T, stringsAsFactors=F, sep=',')
  crispr_minscores <- apply(crispr_depmap_df_all_genes[, 2:ncol(crispr_depmap_df_all_genes)], 1, min)
  names(crispr_minscores) <- crispr_depmap_df_all_genes[, 'gene']
  
  depmap_filename <- paste(depmap_dir, 'D2_combined_gene_dep_scores_all.csv', sep='')
  rnai_depmap_df_all_genes <- read.table(depmap_filename, header=T, stringsAsFactors=F, sep=',')
  rnai_minscores <- apply(rnai_depmap_df_all_genes[, 2:ncol(rnai_depmap_df_all_genes)], 1, min)
  names(rnai_minscores) <- rnai_depmap_df_all_genes[, 'gene']
  
  all_genes <- union(names(rnai_minscores), names(crispr_minscores))
  depmap_min_df_all_genes <- data.frame(gene=all_genes, min_crispr=crispr_minscores[all_genes], min_rnai=rnai_minscores[all_genes])
  depmap_min_df_all_genes[is.na(depmap_min_df_all_genes)] <- 0
  
  return(depmap_min_df_all_genes)
}

# returns a df with 2 cols for INPUT genes: minimum score across all cell lines in RNAi and CRISPR screens
# for efficiency, it takes entire data frames read by calculate_depmap_min_df as input
get_depmap_min_df <- function(input_genes, depmap_min_df_all_genes){
  depmap_min_df <- depmap_min_df_all_genes[match(input_genes, depmap_min_df_all_genes[, 'gene']), c('min_crispr', 'min_rnai')]
  return(depmap_min_df)
}

# returns a vector with avg depmap score values in cancer_type-specific cell lines
# if gene not in depmap list, value is 0
get_depmap_df <- function(input_genes, depmap_dir, depmap_type, cancer_type){
  depmap_filename <- paste(depmap_dir, 'D2_combined_gene_dep_scores_', cancer_type, '.csv', sep='')
  if(depmap_type == 'crispr')
    depmap_filename <- paste(depmap_dir, 'Achilles_gene_effect_', cancer_type, '.csv', sep='')
  
  if(!file.exists(depmap_filename))
    depmap_filename <- gsub(cancer_type, 'all', depmap_filename)
  
  depmap_data <- read.table(depmap_filename, header=T, stringsAsFactors=F, sep=',')
  depmap_data$p_neg_score <- rowSums(depmap_data < 0) / rowSums(depmap_data != 0) # percentage of cell lines in which each gene has -ve value
  
  dp_inds <- match(input_genes, depmap_data$gene)
  
  p_val <- wilcox.test(depmap_data$p_neg_score[dp_inds], depmap_data$p_neg_score[setdiff(1:nrow(depmap_data), dp_inds)], alternative='greater')$p.value # p-value on whether selected genes statistically impact cell survival (i.e. have negative values in a large number of cell lines) more than remaining genes
  print(paste('Mann-Whitney U Test p-value (', depmap_type, '): ', p_val, sep=''))
  
  depmap_df <- depmap_data[dp_inds, c('dep_score', 'p_neg_score')]
  colnames(depmap_df) <- gsub('score', paste('score', depmap_type, sep='_'), colnames(depmap_df))
  
  depmap_df[is.na(depmap_df)] <- 0
  
  return(depmap_df)
}

library(mgsub); library(stats) # replace cancer types to include tissue name
get_cancerMine_vector <- function(input_genes, cancerMine_dir, cancerMine_filename, cancerMine_dictionary_filename, cancer_type){
  cancerMine_data <- read.table(paste(cancerMine_dir, cancerMine_filename, sep=''), header=T, quote="\"", stringsAsFactors=F, fill=T, sep='\t')
  cancerMine_data[, 'cancer_normalized'] <- gsub(' ', '_', cancerMine_data[, 'cancer_normalized'])
  cancerMine_dict <- read.table(paste(cancerMine_dir, cancerMine_dictionary_filename, sep=''), header=T, stringsAsFactors=F, fill=T, sep='\t')
  
  cancerMine_data$cancer_normalized_w_types <- mgsub(cancerMine_data$cancer_normalized, cancerMine_dict$from, cancerMine_dict$to) # add tissue and cancer types (i.e. could be done once and passed as parameter if speedup is needed)
  cancerMine_data <- cancerMine_data[grepl(cancer_type, cancerMine_data$cancer_normalized_w_types), ] # retain data of cancer type only
  cancerMine_data <- aggregate(cancerMine_data[, c('citation_count')], by=list(gene_normalized=cancerMine_data$gene_normalized), sum) # aggregate data frame
  colnames(cancerMine_data) <- c('gene_normalized', 'citation_count')
  
  cancerMine_vector <- cancerMine_data$citation_count[match(input_genes, cancerMine_data$gene_normalized)]
  cancerMine_vector[is.na(cancerMine_vector)] <- 0
  
  return(cancerMine_vector)
}