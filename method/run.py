import ngr
import helper_functions as hp
import evaluation_functions as eval

import numpy as np
from scipy.spatial import distance

import datetime
import uuid
import sys

def main():
    uid = str(uuid.uuid4())[0:8] # unique id
    print('UID: ' +str(uid))

    # get data (matrix and scores)
    initial_score_file = '../base/results/combined_scores_somatic_MC3_feature_files.csv' #'../base/results/test_initial_scores.csv' #
    matrix_file_prefix =  '../ppi/code/results/STRING_converted_matrix' #'../ppi/code/results/test_converted_matrix'
    initial_scores, ppi_matrix, ppi_matrix_index = hp.get_ngr_inputs(initial_score_file, matrix_file_prefix, matrix_normalization_method='none')

    ppi_matrix[ppi_matrix > 0] = 1 # binarize matrix

    # get largest connected component
    ppi_matrix, ppi_matrix_index = hp.get_largest_connected_component(ppi_matrix, ppi_matrix_index)
    print(ppi_matrix)
    print(ppi_matrix_index)
    
    # add new edges
    new_edges_file = '../ppi/curated_pathway_edges/new_test_edges.txt'
    new_edges_suffix = 'sanchez_etal_pathway_edges'
    ppi_matrix = hp.add_new_edges(ppi_matrix, ppi_matrix_index, new_edges_file)
    ppi_matrix = hp.normalize_matrix(ppi_matrix, matrix_file_prefix+'_w_'+new_edges_suffix, normalization='insulated_diffusion') # normalize matrix after adding the new edges
    print(ppi_matrix)

    sys.exit()
    
    # run method
    ngr_scores = ngr.run(Y=initial_scores['score'], W=ppi_matrix, W_index=ppi_matrix_index, initial_scores=initial_scores) # ngr_scores have the order of original initial scores (as passed to this function which matches matrix index order)
    
    # evaluate results
    gold_standard_file = '../gene_lists/code/results/cancerMine_both_depmap_both_values_mean_scores_all_tissue.csv'
    ngr_list, initial_scores = hp.process_ngr_results(initial_scores, ngr_scores, uid=uid) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores
    eval.evaluate_results((initial_scores, ngr_list), gold_standard_file, ('Initial', 'NGR'), comp_measures=('AO', 'weighted_tau', 'manhattan_dist', 'tau', 'RBO', 'euclidean_dist'))

main()
