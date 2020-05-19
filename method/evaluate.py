import helper_functions as hp
import evaluation_functions as eval

import numpy as np
import pickle

## evaluate results
# A. Whole PPI
variant_type = 'somatic_MC3'; cancer_type = 'BRCA'; ppi_network='HumanNetv2';
gold_standard_file = '../gene_lists/code/results/cancerMine_both_depmap_both_values_mean_scores_breast_tissue.csv'

# load whole-ppi results
uid = 'ab420746'; output_dir='results/lcc_whole_network/';
output_pickle = output_dir+variant_type+'_'+cancer_type+'_'+ppi_network+'_lcc_'+uid+'_score_matrix.pkl'
ngr_score_matrix = pickle.load(open(output_pickle, "rb"))

# load whole ppi;
ppi_matrix_prefix = '../ppi/code/results/'+ppi_network+'_converted_matrix'
genomics_matrix_prefix = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'+variant_type+'_'+cancer_type+'_matrix'

# read whole-ppi index; IMPORTANT SET VALUE OF PARAMETER ppi_cancer_type TO '' TO GET INDEX OF WHOLE PPI FOR RESULT ASSESSMENT
ppi_matrix, ppi_matrix_index, extended_ppi_matrix_prefix, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = hp.get_ngr_inputs(genomics_matrix_prefix, ppi_matrix_prefix, genomics_matrix_transformed=True, largest_cc=True, ppi_cancer_type='', matrix_normalization_method='insulated_diffusion')

print(ngr_score_matrix.shape)
print(ppi_matrix_index.shape)
print('--------------------')

print(genomics_matrix_row_index.shape)
print(genomics_matrix_col_index.shape)
print(genomics_matrix.shape)
print('--------------------')

print(extended_ppi_matrix_prefix)

#B. Tissue-specific PPI

# load tissue-specific-ppi results
uid2 = '55b78775'; output_dir2='results/';
output_pickle2 = output_dir2+variant_type+'_'+cancer_type+'_'+ppi_network+'_lcc_'+cancer_type.lower()+'_'+uid2+'_score_matrix.pkl'
ngr_score_matrix2 = pickle.load(open(output_pickle2, "rb"))

# load tissue-specific ppi;
ppi_matrix_prefix2 = '../ppi/code/results/'+ppi_network+'_converted_matrix'
genomics_matrix_prefix2 = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'+variant_type+'_'+cancer_type+'_matrix'

# read tissue-specific-ppi index; IMPORTANT SET VALUE OF PARAMETER ppi_cancer_type TO cancer_type TO FILTER TISSUE-SPECIFC UNEXPRESSED GENES OUT FOR RESULT ASSESSMENT
ppi_matrix2, ppi_matrix_index2, extended_ppi_matrix_prefix2, genomics_matrix2, genomics_matrix_row_index2, genomics_matrix_col_index2 = hp.get_ngr_inputs(genomics_matrix_prefix2, ppi_matrix_prefix2, genomics_matrix_transformed=True, largest_cc=True, ppi_cancer_type=cancer_type, matrix_normalization_method='insulated_diffusion')

print(ngr_score_matrix2.shape)
print(ppi_matrix_index2.shape)
print('--------------------')

print(genomics_matrix_row_index2.shape)
print(genomics_matrix_col_index2.shape)
print(genomics_matrix2.shape)
print('--------------------')

print(extended_ppi_matrix_prefix2)

# C. Load scores and evaluate
ngr_scores = hp.process_ngr_results(ngr_score_matrix, ppi_matrix_index, uid=uid, output_dir=output_dir) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores
ngr_scores2 = hp.process_ngr_results(ngr_score_matrix2, ppi_matrix_index2, uid=uid2, output_dir=output_dir2) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores

eval.evaluate_results((ngr_scores2, ngr_scores), gold_standard_file, ('Tissue-PPI', 'Whole-PPI'), cancer_type=cancer_type, variant_type=variant_type, comp_measures=('Accuracy', 'F1', 'Precision'))