source('plotting.R')

library(dplyr)
library(scales)
library(ggplot2)

variant_type <- 'somatic_MC3'
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0 && args[1] == '-v')
    variant_type <- args[2]

process_columns <- function(x, gene_expression_col=F){
  xnz <- x[x != 0] # non zero values
  xnz <- abs(xnz) # absolute values; applicable to all columns but more specific to gene expression columns (over- = under-expression, 0 = no diff exp.)
  
  if(gene_expression_col == F)
    xnz <- log(xnz)
  
  xnz <- rescale(xnz) + 0.0000001 # minmax normalization (rescale of scales package) + epsilon
  x[x!=0] <- xnz # set non-zero values back in original column
  return(x)
}

feature_files_dir <- 'feature_file_lists/'
feature_files_basename <- paste(variant_type, 'feature', 'files', sep='_')
feature_files_filename <- paste(feature_files_dir, feature_files_basename, '.txt', sep='')
feature_files <- read.table(feature_files_filename, stringsAsFactors=F)[, 1] # read the list of files

# All feature files have: gene names in 1st column, feature value of interest in 2nd
feature_data_combined <- read.table(feature_files[1], header=T, sep=',', stringsAsFactors=F)[, 1:2]
print(paste('Reading', feature_files[1], 'Done.', sep=' '))

for(i in 2:length(feature_files)){
  feature_data <- read.table(feature_files[i], header=T, sep=',', stringsAsFactors=F)
  feature_data <- feature_data[which(!duplicated(feature_data)), ] # remove duplicates to be on the safe side

  if(nrow(feature_data) > 0){
    feature_data_combined <- full_join(feature_data_combined, feature_data[, 1:2], by="gene") # merge files
    print(paste('Processing', feature_files[i], 'Done.', sep=' '))
  } else{
    feature_data_combined[, ncol(feature_data_combined)+1] <- rep(1, nrow(feature_data_combined))
    print(paste('Empty file: ', feature_files[i], '| Added an dummy column of 1s. Moving on.', sep=' '))
  }
}

print(head(feature_data_combined))
print(dim(feature_data_combined))

feature_data_combined[is.na(feature_data_combined)] <- 0 # replace NAs with 0's
colnames(feature_data_combined)[c(2,3,6,7,8)] <- c('varpergene_normalized', 'samplepergene_normalized', 'nonsynonym_normalized', 'stopgain_normalized', 'frameshift_normalized')

non_gene_expression_col_start <- 2; non_gene_expression_col_end <- 13
feature_data_combined[, non_gene_expression_col_start:non_gene_expression_col_end] <- apply(feature_data_combined[, non_gene_expression_col_start:non_gene_expression_col_end], 2, process_columns) # process non-gene expression columns

gene_expression_col_start <- 14; gene_expression_col_end <- 31
feature_data_combined[, gene_expression_col_start:gene_expression_col_end] <- apply(feature_data_combined[, gene_expression_col_start:gene_expression_col_end], 2, process_columns, gene_expression_col=T) # process gene expression columns

tissue <- 'all'
if(tissue == 'all')
  feature_data_combined$logFC_score <- rescale(apply(feature_data_combined[, gene_expression_col_start:gene_expression_col_end], 1, mean)) # average expression value
# else statement for tissue-specific aggregation

# known gene list membership scaled from 1 to to 0.2
known_gene_list_features <- c(11:13)
feature_data_combined[, known_gene_list_features] <- (feature_data_combined[, known_gene_list_features]/5)

# select features
selected_features <- c(1:7, 9, 10, 32) # exclude frameshift, known lists, and expression raw values for now (aggregate expression value at index 32 is used)
feature_data_combined <- feature_data_combined[, selected_features]
  
feature_data_combined$combined_score <- apply(feature_data_combined[2:ncol(feature_data_combined)], 1, mean) # calculate combined scores
feature_data_combined <- feature_data_combined[order(feature_data_combined$combined_score, decreasing=TRUE), ] # sort by combined score

print('Top genes:')
head(feature_data_combined[, c('gene', 'combined_score')], 10)

print('Bottom genes:')
tail(feature_data_combined[, c('gene', 'combined_score')], 10)

output_dir <- 'results/'
#save(feature_data_combined, file=paste(output_dir, 'combined_scores_', feature_files_basename, '.RData', sep=''))
write.table(feature_data_combined[, c('gene', 'combined_score')], file=paste(output_dir, 'combined_scores_', feature_files_basename, '.csv', sep=''), col.names=T, row.names=F, sep=',', quote=F)

# Visualization
if(length(args) > 2){
	if(args[3] == '-p')
		generate_histogram_from_list(feature_data_combined$combined_score, x_label='NGR Combined Score', title='', bins=10)
}


# Evaluation
evaluate_list <- function(param_list, target_list, evaluation_type='top_percent', top_percent=0.1){
  if(evaluation_type == 'top_percent'){ # check the ration of top_percent of target list are in the top part of param_list
    n_top_elements <- min(floor(top_percent * length(target_list)), length(param_list))
    print(n_top_elements)
    n_intersection <- length(intersect(param_list[1:n_top_elements], target_list[1:n_top_elements]))
    print(n_intersection)
    intersection_ratio <- (n_intersection / n_top_elements)
    print(paste('Intersection ratio: ', intersection_ratio, ' (', n_intersection, '/', n_top_elements, ').', sep=''))
    return(intersection_ratio)
  }
}

#target_file <- '../gene_lists/dependency/D2_combined_gene_dep_scores_processed.csv'
#target_data <- read.table(target_file, sep=',', header=T, stringsAsFactors=F)
#target_data <- target_data[order(abs(target_data[, 'dep_score']), decreasing=T), ]
#evaluate_list(param_list=feature_data_combined[, 'gene'], target_list=target_data[, 'gene'])
