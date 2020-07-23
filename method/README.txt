# to generate whole mobility lists, create results/mobility_lists/whole_lists dir once then locally run the below with " Mobility List Batch Generation" uncommented in evaluate.py:
python evaluate.py > results/mobility_lists/whole_lists/mobility_lists.log

# to generate NGR matrix result on 18 cancer types, in method/, run:
sbatch run_ngr_all_cancer_types.sh

# to generate complete mobility lists on 18 cancer types, in method/, run:
sbatch generate_all_mobility_lists.sh

# to generate prioritized mobility lists used in results, in ../analysis/:
Run mobility_lists_analysis.Rmd
+ capture output manually and save it in NGR/method/results/mobility_lists/output_w_gene_list_details.log; the output includes known genes captured by upward mobility and number of upward mobility genes across lists.

# to generate hypergemetric results, uncomment the hypergemetric results block in evalaute.py then in method/, run:
python evaluate.py > results/hypergemetric_test_results.txt