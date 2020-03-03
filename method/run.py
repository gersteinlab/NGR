import ngr
import helper_functions as hp
import evaluation_functions as eval

import numpy as np
from scipy.spatial import distance

import datetime
import uuid

def main():
    uid = str(uuid.uuid4())[0:8] # unique id
    print('UID: ' +str(uid))

    # A. Method Execution
    initial_score_file = '../base/results/combined_scores_somatic_MC3_feature_files.csv'
    matrix_file_prefix = '../ppi/code/results/STRING_converted_matrix' #test_matrix2
    initial_scores, ppi_matrix, ppi_matrix_index = hp.get_ngr_inputs(initial_score_file, matrix_file_prefix)

    print(initial_scores[1:5])
    print(initial_scores.shape)
    
    ngr_scores = ngr.run(Y=initial_scores['score'], W=ppi_matrix, W_index=ppi_matrix_index, initial_scores=initial_scores) # ngr_scores have the order of original initial scores (as passed to this function which matches matrix index order)
    
    # B. Evaluation
    gold_standard_file = '../gene_lists/code/results/cancerMine_both_depmap_both_values_mean_scores_all_tissue.csv'
    ngr_list, initial_scores = hp.process_ngr_results(initial_scores, ngr_scores, uid=uid) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores
    eval.evaluate_results((initial_scores, ngr_list), gold_standard_file, ('Initial', 'NGR'), comp_measures=('AO', 'weighted_tau', 'manhattan_dist', 'tau', 'RBO', 'euclidean_dist'))

main()
