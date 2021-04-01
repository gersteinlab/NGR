import helper_functions as hp
import evaluation_functions as eval

import numpy as np
import datetime

# If revisited: 
# Might want to handle ties in terms of same initial score: i.e. assign same initial score the same rank and adjust T beta value accordingly
# in pushed code, initial and ngr/umg ranks start at 1 (hence +1 to initial_ranks and postprop_ranks in calculate_mobiltiy()), in current code used to generate paper results it starts at 0 (hence no +1 in generate_gene_mobility_lists())

# Y is the vector of initial scores
# W is the weight matrix, i.e. PPI network
# W_index is the ordered list of genes as they appear in cols and rows
def run(M, W, alpha=0.8):
    print('\nNGR diffusion starts...\n')
    i = 0; epsilon = 1 / np.power(10, 6); max_iterations = 350
    M0 = M; M_delta = epsilon

    print('{0}\nPropagation in progress...'.format(datetime.datetime.now()))
    while M_delta >= epsilon and i < max_iterations:
        previous = M

        M = (alpha * np.dot(W, M)) + ((1 - alpha) * M0) # propagation step      
        M_delta = np.linalg.norm(M - previous)
        print('\tIter. {0} complete. Heat difference: {1}'.format(i+1, M_delta))
        
        i += 1

    print('Propagation done: {0} iterations.'.format(i))
    print(datetime.datetime.now())

    return M
