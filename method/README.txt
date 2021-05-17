# to generate whole mobility lists, create results/mobility_lists/whole_lists dir once then locally run the below with " Mobility List Batch Generation" uncommented in evaluate.py:
python evaluate.py > results/mobility_lists/whole_lists/mobility_lists.log

# to generate NGR matrix result on 17 cancer types, in method/, run:
sbatch run_ngr_all_cancer_types.sh

# to generate complete mobility lists on 17 cancer types, in method/, run:
sbatch generate_all_mobility_lists.sh

# to generate prioritized mobility lists used in results, in ../analysis/:
Run mobility_lists_analysis.Rmd

# to generate Mann-Whitney U one-sided test results for ranks, uncomment the generate_batch_pvalues() block in evalaute.py then in method/, run:
python evaluate.py > results/mann_whitney_U_onesided_test_of_ranks_results.txt
