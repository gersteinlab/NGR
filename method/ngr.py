import helper_functions as hp
import evaluation_functions as eval

import numpy as np
import datetime

# Y is the vector of initial scores
# W is the weight matrix, i.e. PPI network
# W_index is the ordered list of genes as they appear in cols and rows
def run(Y, W, W_index, initial_scores=""):
    Y = hp.totality_scale(Y.reshape(Y.shape[0], 1)) # scale values so they all sum to 1    
    S = np.ndarray.copy(Y)
    alpha = 0.8

    print('Initial scores:{}'.format(initial_scores))
    print('Y: {}'.format(Y.reshape(Y.shape[0], )))

    print('\nNGR diffusion starts...\n')
    n_iterations = 30; i = 0; S_total_heat = 0
    while i < n_iterations:
        #print('Scores at beg of iteration {0}: {1}'.format(i, S.T))
                
        Sm = np.repeat(S, S.shape[0], axis=1)
        Sm = (W > 0) * Sm # Step 1: W-based ReLU on S                
        Sm = np.apply_along_axis(hp.nonzero_softmax, axis=0, arr=Sm, inverted_x=True) # Step 2: row-wise softmax
        
        #print(Sm)
        
        S = ((1 - alpha) * Y) + (alpha * np.matmul(Sm, S)) # Step 3: S = (1-a)*(Y) + (a)*Sm.S
        
        #S = hp.totality_scale(S) # re-normalize after each step; safety step 

        S_total_heat_prev = S_total_heat
        S_total_heat = np.sum(S)
        heat_diff = np.abs(S_total_heat - S_total_heat_prev)

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

    print('Final scores:{}'.format(S.reshape(S.shape[0], )))

    return S.reshape(S.shape[0], )
