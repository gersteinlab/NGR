import ngr
import helper_functions as hp
import evaluation_functions as eval

import numpy as np
from scipy.spatial import distance

import datetime
import uuid
import pickle
import sys

def main():
    uid = str(uuid.uuid4())[0:8] # unique id
    print('UID: ' +str(uid))

    variant_type = 'somatic_MC3'; cancer_type = 'BRCA'; ppi_network='STRING'; output_dir='results/'; output_filename_suffix = 'lcc'    
    ppi_matrix_prefix = '../ppi/code/results/'+ppi_network+'_converted_matrix'
    genomics_matrix_prefix = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'+variant_type+'_'+cancer_type+'_matrix'

      
    # get data (ppi_matrix and scores)
    ppi_matrix, ppi_matrix_index, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = hp.get_ngr_inputs(genomics_matrix_prefix, ppi_matrix_prefix, largest_cc=True, matrix_normalization_method='insulated_diffusion')
    print('PPI matrix: {0}, Input matrix: {1}.\n'.format(ppi_matrix.shape, genomics_matrix.shape))
    '''
    genomics_matrix[genomics_matrix > 0] = 1 # binarize genomics matrix
    
    # add new edges and normalize matrix
    #ppi_matrix[ppi_matrix > 0] = 1 # binarize matrix before adding edges
    #new_edges_file = '../ppi/curated_pathway_edges/SanchezVega-etal-pathway_edges.txt'
    #new_edges_suffix = 'sanchez_etal_pathway_edges'
    #output_filename_suffix += '_w_'+new_edges_suffix
     
    #ppi_matrix = hp.add_new_edges(ppi_matrix, ppi_matrix_index, new_edges_file)
    #ppi_matrix = hp.normalize_matrix(ppi_matrix, ppi_matrix_prefix+'_'+output_filename_suffix, normalization='insulated_diffusion') # normalize matrix after adding new edges
    
    # run method
    ngr_score_matrix = ngr.run(X=genomics_matrix, W=ppi_matrix)
    pickle.dump(ngr_score_matrix, open('_'.join([output_dir+variant_type, cancer_type, ppi_network, output_filename_suffix, uid, 'score_matrix.pkl']), 'wb'))
    '''
    
    output_pickle = output_dir+'somatic_MC3_BRCA_HumanNetv2_lcc_ab420746_score_matrix.pkl'
    uid = 'ab420746'
    ngr_score_matrix = pickle.load(open(output_pickle, "rb"))

    # evaluate results
    gold_standard_file = '../gene_lists/code/results/cancerMine_both_depmap_both_values_mean_scores_breast_tissue.csv'
    ngr_scores = hp.process_ngr_results(ngr_score_matrix, ppi_matrix_index, uid=uid, output_dir=output_dir) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores
    
    eval.evaluate_results((ngr_scores, ngr_scores), gold_standard_file, ('NGR', 'NGR'), N=1500, comp_measures=('AUC', 'Hypergeomtric_test'))
    sys.exit()
    
main()
