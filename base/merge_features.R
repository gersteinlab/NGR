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

library(dplyr)

cancer_type = "" # default is empty string and implies pancancer analysis
if(cancer_type == ""){
  feature_files_filename <- paste("feature_file_lists/", cancer_type, "feature_files.txt", sep="")
} else
  feature_files_filename <- paste("feature_file_lists/", cancer_type, "_feature_files.txt", sep="")

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