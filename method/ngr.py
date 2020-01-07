import helper_functions as hp

import numpy as np
import datetime

# Y is the vector of initial scores
# W is the weight matrix, i.e. PPI network
# W_index is the ordered list of genes as they appear in cols and rows
def run(Y, W, W_index):
    Y = hp.totality_scale(Y.reshape(Y.shape[0], 1)) # scale values so they all sum to 1    
    S = np.ndarray.copy(Y)
    alpha = 0.8

    print('\nNGR diffusion starts...\n')    
    n_iterations = 2; i = 0; S_total_heat = 0
    while i < n_iterations:
        Sm = np.repeat(S, S.shape[0], axis=1)
        Sm = (W > 0) * Sm # Step 1: W-based ReLU on S 
        Sm = np.apply_along_axis(hp.nonzero_softmax, axis=1, arr=Sm) # Step 2: row-wise softmax
        S = ((1 - alpha) * Y) + np.matmul(Sm, S) # Step 3: S = (1-a)*(Y) + (a)*Sm.S
        
        print(S)
        S_total_heat_prev = S_total_heat
        S_total_heat = np.sum(S)
        heat_diff = np.abs(S_total_heat - S_total_heat_prev)
        print('Iteration {0} complete. Heat difference: {1}'.format(i, heat_diff))

        i += 1

    # start here: in-place, 0-based softmax,
    print('\nNGR diffusion done.')
    print(datetime.datetime.now())

    return S.reshape(S.shape[0], )
