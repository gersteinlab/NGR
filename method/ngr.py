import helper_functions as hp
import numpy as np

def main():
    initial_score_file = '../base/results/combined_scores_somatic_MC3_feature_files.csv'
    initial_scores = hp.read_scores(initial_score_file, sorted=True)
    initial_scores = initial_scores[np.argsort(initial_scores['gene'])]
    print('Initial scores read successfully.')

    # numpy array of voids
    ppi_matrix_file_prefix = '../ppi/code/results/test_matrix2'
    ppi_matrix, ppi_matrix_index = hp.read_matrix(ppi_matrix_file_prefix, sorted=True, selected_genes=initial_scores['gene'])
    print('PPI matrix read successfully.')
    
    # reorder initial scores to match order of (NxN) matrix cols and rows 
    initial_scores_new_order = np.searchsorted(initial_scores['gene'], ppi_matrix_index) # select only common genes for initial_scores
    initial_scores = initial_scores[initial_scores_new_order]

    # weight matrix normalization
    matrix_normalization_method = ''
    matrix_scaling = True
    
    ppi_matrix = hp.normalize_matrix(ppi_matrix, scaling=matrix_scaling, normalization=matrix_normalization_method)
    print(ppi_matrix)
    print('Weight matrix normalized successfully.')
    
    # callling NGR
    scores = hp.ngr(Y=initial_scores['score'], W=ppi_matrix, W_index=ppi_matrix_index)
    print(scores)
    print(initial_scores)

main()
