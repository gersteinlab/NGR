# library(argparser)
# 
# parser <- arg_parser(description="Parsing command argument to aggregate features.")
# 
# # specify our desired options 
# # by default ArgumentParser will add an help option 
# add_argument(parser, "--cancer_type", default="",  help="Cancer type for which features are parsed.")
# argv <- parse_args(parser)
# 
# print(argv)

source('plotting.R')

library(dplyr)
library(scales)
library(argparser)

parser <- arg_parser('argument parser')
parser <- add_argument(parser, "-v", default="somatic_MC3", help="Variation type: somatic_MC3 or germline")
parser <- add_argument(parser, "-p", default=FALSE, help="Plotting flag")
args = parse_args(parser)

variant_type <- "somatic_MC3"
if(args$v != "")
  variant_type <- args$v

feature_files_dir <- "feature_file_lists/"
feature_files_basename <- paste(variant_type, "feature", "files", sep="_")
feature_files_filename <- paste(feature_files_dir, feature_files_basename, ".txt", sep="")
feature_files <- as.character(read.table(feature_files_filename)[, 1]) # read the list of files

# All feature files have: gene names in 1st column, feature value of interest in 2nd
feature_data_combined <- read.table(feature_files[1], header=T, sep=",")[, 1:2]
print(paste("Reading", feature_files[1], "Done.", sep=" "))

for(i in 2:length(feature_files)){
  feature_data <- read.table(feature_files[i], header=T, sep=",")
  feature_data_combined <- full_join(feature_data_combined, feature_data[, 1:2]) # merge files
  print(paste("Processing", feature_files[i], "Done.", sep=" "))
}

print(head(feature_data_combined))
print(dim(feature_data_combined))

feature_data_combined[is.na(feature_data_combined)] <- 0 # replace NAs with 0's
feature_data_combined[, 2:ncol(feature_data_combined)] <- apply(feature_data_combined[, 2:ncol(feature_data_combined)], 2, rescale) # min-max scaling per column
feature_data_combined$combined_score <- apply(feature_data_combined[2:ncol(feature_data_combined)], 1, mean) # calculate combined scores
feature_data_combined <- feature_data_combined[order(feature_data_combined$combined_score, decreasing=TRUE), ] # sort by combined score

print("Top genes:")
head(feature_data_combined[, c("gene", "combined_score")], 10)

print("Bottom genes:")
tail(feature_data_combined[, c("gene", "combined_score")], 10)

save(feature_data_combined, file=paste("combined_scores_", feature_files_basename, ".RData", sep=""))

if(args$p)
  generate_histogram_from_list(feature_data_combined$combined_score, x_label='NGR Combined Score', title='', bins=10)
