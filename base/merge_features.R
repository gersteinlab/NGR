# library(argparser)
# 
# parser <- arg_parser(description='Parsing command argument to aggregate features.')
# 
# # specify our desired options 
# # by default ArgumentParser will add an help option 
# add_argument(parser, '--cancer_type', default='',  help='Cancer type for which features are parsed.')
# argv <- parse_args(parser)
# 
# print(argv)

source('plotting.R')

library(dplyr)
library(scales)
library(argparser)
library(ggplot2)

parser <- arg_parser('argument parser')
parser <- add_argument(parser, '-v', default='somatic_MC3', help='Variation type: somatic_MC3 or germline')
parser <- add_argument(parser, '-p', default=FALSE, help='Plotting flag')
args = parse_args(parser)

variant_type <- args$v

feature_files_dir <- 'feature_file_lists/'
feature_files_basename <- paste(variant_type, 'feature', 'files', sep='_')
feature_files_filename <- paste(feature_files_dir, feature_files_basename, '.txt', sep='')
feature_files <- as.character(read.table(feature_files_filename)[, 1]) # read the list of files

# All feature files have: gene names in 1st column, feature value of interest in 2nd
feature_data_combined <- read.table(feature_files[1], header=T, sep=',')[, 1:2]
print(paste('Reading', feature_files[1], 'Done.', sep=' '))

for(i in 2:length(feature_files)){
  feature_data <- read.table(feature_files[i], header=T, sep=',')
  if(nrow(feature_data) > 0){
    feature_data_combined <- full_join(feature_data_combined, feature_data[, 1:2], by="gene") # merge files
    print(paste('Processing', feature_files[i], 'Done.', sep=' '))
  } else{
    print(paste('Empty file: ', feature_files[i], '| Moving on.', sep=' '))
  }
}

print(head(feature_data_combined))
print(dim(feature_data_combined))

feature_data_combined[is.na(feature_data_combined)] <- 0 # replace NAs with 0's
feature_data_combined[, 2:ncol(feature_data_combined)] <- apply(feature_data_combined[, 2:ncol(feature_data_combined)], 2, rescale) # min-max scaling per column
feature_data_combined$combined_score <- apply(feature_data_combined[2:ncol(feature_data_combined)], 1, mean) # calculate combined scores
feature_data_combined <- feature_data_combined[order(feature_data_combined$combined_score, decreasing=TRUE), ] # sort by combined score

print('Top genes:')
head(feature_data_combined[, c('gene', 'combined_score')], 10)

print('Bottom genes:')
tail(feature_data_combined[, c('gene', 'combined_score')], 10)

output_dir <- 'results/'
save(feature_data_combined, file=paste(output_dir, 'combined_scores_', feature_files_basename, '.RData', sep=''))
write.table(feature_data_combined[, c('gene', 'combined_score')], file=paste(output_dir, 'combined_scores_', feature_files_basename, '.csv', sep=''), col.names=T, row.names=F, sep=',', quote=F)

# Visualization
if(args$p)
  generate_histogram_from_list(feature_data_combined$combined_score, x_label='NGR Combined Score', title='', bins=10)

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

target_file <- '../gene_lists/dependency/D2_combined_gene_dep_scores_processed.csv'
target_data <- read.table(target_file, sep=',', header=T)
target_data <- target_data[order(abs(target_data[, 'dep_score']), decreasing=T), ]
evaluate_list(param_list=feature_data_combined[, 'gene'], target_list=target_data[, 'gene'])

