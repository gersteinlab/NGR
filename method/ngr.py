import helper_functions as hp
import evaluation_functions as eval

import numpy as np
import datetime

# Y is the vector of initial scores
# W is the weight matrix, i.e. PPI network
# W_index is the ordered list of genes as they appear in cols and rows
def run(X, W, alpha=0.8):
    print('\nNGR diffusion starts...\n')
    n_iterations = 450; i = 0; 
    X0 = X; X_total_heat = 0; heat_diff=0.00001
    
    while i < n_iterations and heat_diff > 0:
        X = (alpha * np.dot(W, X)) + ((1 - alpha) * X0)
        
        #Sm = (W > 0) * Sm # Step 1: W-based ReLU on S                
        #Sm = np.apply_along_axis(hp.nonzero_softmax, axis=0, arr=Sm, inverted_x=True) # Step 2: row-wise softmax
        #S = hp.totality_scale(S) # re-normalize after each step; safety step 

        X_total_heat_prev = X_total_heat
        X_total_heat = np.sum(X)
        heat_diff = np.abs(X_total_heat - X_total_heat_prev)

        #print('Scores at end of iteration {0}: {1}'.format(i, S.T))
        
        # To remove this final block and the initial_scores parameter in function header in the final code
        #gold_standard_file = '../gene_lists/code/results/cancerMine_both_depmap_both_values_mean_scores_all_tissue.csv'
        #ngr_list, initial_scores = hp.process_ngr_results(initial_scores, S.reshape(S.shape[0], ), save_files=False) # ngr and initial scores in both structured arrays (i.e. called lists here) are re-ordered w.r.t. to their respective decreasing order of scores
        #eval.evaluate_results((initial_scores, ngr_list), gold_standard_file, ('Initial', 'NGR'),  N=3000, comp_measures=('AO', 'weighted_tau', 'manhattan_dist', 'tau', 'RBO', 'euclidean_dist'))

        print('\tIter. {0} complete. Heat difference: {1}'.format(i+1, heat_diff))
        i += 1

    # start here: in-place, 0-based softmax,
    print('\nNGR diffusion done.')
    print(datetime.datetime.now())

    print(X)
    
    return X
