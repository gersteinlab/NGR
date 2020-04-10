import numpy as np
from scipy.special import softmax
import networkx as nx 

import datetime
import pickle
import uuid

from os import path

epsilon = 0.0000001

## Helper functions to method execution

# calculate softmax of non-zero values
def nonzero_softmax(x, inverted_x=False):
    x_temp = x
    nonzero_indices = np.asarray(np.where(x_temp != 0))

    if inverted_x: # invert non-zero values (nzx = totality_scale(1-nzx) so that lower values get higher scores yet sum of vector is still 1)
        x_temp[nonzero_indices] = totality_scale(1-x_temp[nonzero_indices])

    if nonzero_indices.size > 0:
        softmax_values = softmax(x_temp[nonzero_indices])
        x_temp[nonzero_indices] = softmax_values
        
    return x_temp

# scaling a vector so its values sum to 1
def totality_scale(values):
    total = values.sum() + epsilon
    return (values/total)

# score files have two columns: gene names and scores
def read_scores(file_name, sorted=True):
    scores = np.loadtxt(file_name, delimiter=',', skiprows=1, 
                        dtype= {'names': ('gene', 'score'), 'formats': ('U20', 'f')})
    
    if sorted:
        sorted_order = np.argsort(scores['gene'])
        scores = scores[sorted_order]

    print('Scores read successfully.')
    print(datetime.datetime.now())
    
    return scores

# reads the PPI network matrix
def read_matrix(file_name_prefix, sorted=True, selected_genes=[]):
    
    matrix_pickle_file_name = file_name_prefix+'.pkl'
    matrix_pickle_index_file_name = file_name_prefix+'_index.pkl'

    if path.exists(matrix_pickle_file_name) and path.exists(matrix_pickle_index_file_name):
        matrix = pickle.load(open(matrix_pickle_file_name, 'rb'))
        matrix_index = pickle.load(open(matrix_pickle_index_file_name, 'rb'))
        print('\nMatrix and index reading done. Existing pickles loaded.')
    else:
        matrix_file_name = file_name_prefix+'.txt'
        matrix_index_file_name = file_name_prefix+'_index.txt'
        
        matrix = np.loadtxt(matrix_file_name, dtype='f') # pickle matrix
        with open(matrix_pickle_file_name, 'wb') as f:
            pickle.dump(matrix, f)

        matrix_index = np.loadtxt(matrix_index_file_name, dtype='U25') # pickle matrix index
        #with open(matrix_pickle_index_file_name, 'wb') as f:
        #    pickle.dump(matrix_index, f)

        print('\nMatrix and index reading done. Pickles saved.')
    
    print(datetime.datetime.now())

    if sorted:
        sorted_order = np.argsort(matrix_index)
        matrix = matrix[sorted_order, :]
        matrix = matrix[:, sorted_order]
        print('Matrix row-col sorting done.')
 
    if (len(selected_genes) > 0): # select genes in rows and columns
        #if(np.diff(selected_genes) >= 0 or np.diff(matrix_index) >= 0):
        #    raise Exception('Warning: selected_genes and/or matrix index are not sorted.')

        common_genes = np.intersect1d(matrix_index, selected_genes)
        common_genes_in_matrix_index = np.searchsorted(matrix_index, common_genes) # indices of common genes in the matrix

        #print(selected_genes)
        #print(common_genes_in_matrix_index)
        
        matrix = matrix[:, common_genes_in_matrix_index] # reorder matrix rows and cols
        matrix = matrix[common_genes_in_matrix_index, :]
        
        matrix_index = matrix_index[common_genes_in_matrix_index] # reorder index

    print('\nMatrix and index selection done.')
    print(datetime.datetime.now())
    print('Final matrix dimensions: ' + str(matrix.shape))
    
    return matrix, matrix_index

# returns an updated matrix including new edges from input edge list 
def add_new_edges(matrix, matrix_index, input_edgelist_file, edge_type='directed'):
    updated_matrix = matrix

    # convert index to a dictionary for fast processing
    dmatrix_index = {};
    for i in range(len(matrix_index)):
        dmatrix_index[matrix_index[i]] = i
        
    # read and add edges to matrix
    separator = " "
    with open(input_edgelist_file, "r") as f:
        for line in f:
            if line.strip() != "" and line[0] != "#": # ignore empty and comment lines
                edge = line.strip().split(separator)
                
                if edge[0] in dmatrix_index and edge[1] in dmatrix_index: # check if both nodes are in the graph
                    source_ind = dmatrix_index[edge[0]]
                    dest_ind = dmatrix_index[edge[1]]
                    edge_weight = abs(np.float(edge[2]))
                    
                    updated_matrix[dest_ind][source_ind] += edge_weight # col-to-row matrix graph representation
    
                    if edge_type == 'undirected':
                        updated_matrix[source_ind][dest_ind] += edge_weight
    
    return updated_matrix
    
def get_largest_connected_component(matrix, matrix_index):
    G = nx.from_numpy_matrix(matrix)
    CCs_obj = nx.connected_components(G)
    
    CCs = []; CC_len = []
    for c in CCs_obj:
        CCs.append(np.array(list(c)))
        CC_len.append(len(c))
    
    largest_CC = CCs[np.argmax(CC_len)]

    print('Number of nodes in largest connected component: {0}'.format(len(largest_CC)))
    print('Indices in largest connected component: {0}'.format(largest_CC))
    print('Nodes in largest connected component: {0}'.format(matrix_index[largest_CC]))
    
    matrix = matrix[largest_CC, :]
    matrix = matrix[:, largest_CC]
    matrix_index = matrix_index[largest_CC]
    
    return matrix, matrix_index # return largest connected component (indices of its vertices)
    
def get_ngr_inputs(initial_score_file, matrix_file_prefix, matrix_normalization_method='insulated_diffusion'):
    initial_scores = read_scores(initial_score_file, sorted=True) # numpy structured array
    
    ppi_matrix, ppi_matrix_index = read_matrix(matrix_file_prefix, sorted=True)

    # reorder initial scores to match order of (NxN) matrix cols and rows; reordered list of scores to be provided to the method to act as Y (i.e. scores) vector
    initial_scores = initial_scores[np.argsort(initial_scores['gene'])] # this step is needed for searchsorted() to work
    initial_scores_new_order = np.searchsorted(initial_scores['gene'], ppi_matrix_index) # select only common genes for initial_scores

    # start here
    # subtract matrix index from list of genes with initial scores
    # remove genes in initial scores list not in ppi
    # search for initial scores genes in updated list and assign them their score according to their order in matrix index; all other genes assigned a score of zero
    
    print(initial_scores_new_order)
    initial_scores = initial_scores[initial_scores_new_order]
    #print(initial_scores['gene'])

    # filter out genes without scores in the ppi matrix
    if matrix_normalization_method != 'none':
        ppi_matrix = normalize_matrix(ppi_matrix, matrix_file_prefix, normalization=matrix_normalization_method) # normalize weight matrix
    
    return initial_scores, ppi_matrix, ppi_matrix_index

# initial_scores is a structured array with two columns: 'gene' and 'score'
# ngr_scores is a numpy arrays of scores
# processes the output of NGR and returns two structured arrays to facilitate evaluation; each array has two columns, 'gene' and 'score,' and its rows are sorted in decreasing order of score
def process_ngr_results(initial_scores, ngr_scores, save_files=True, uid=''):
    # sort result scores in decreasing order 
    ngr_order = np.flip(np.argsort(ngr_scores)) # indices of genes sorted by decreasing order of diffused scores
    ngr_genes = initial_scores['gene'][ngr_order].flatten() # ngr scores are updated ones per order of genes in initial_scores['gene']
    ngr_scores = ngr_scores[ngr_order].flatten()
    ngr_list = np.array(list(zip(ngr_genes, ngr_scores)), dtype=[('gene', 'U25'), ('score', 'f8')]) # create structured array for ngr genes and scores
    
    # sort initial scores in decreasing order 
    initial_scores_order = np.flip(np.argsort(initial_scores['score']))
    initial_scores_sorted = initial_scores['score'][initial_scores_order]
    initial_genes_sorted = initial_scores['gene'][initial_scores_order]
    
    if save_files:
        if uid != '':
            uid = str(uid)+'_'
            
        np.savetxt(uid+'ngr_results.csv', np.stack((ngr_genes, ngr_scores), axis=1), fmt='%s', delimiter=',')
        np.savetxt(uid+'initial_results.csv', np.stack((initial_genes_sorted, initial_scores_sorted), axis=1), fmt='%s', delimiter=',') # initial scores sorted in decreasing order of score value
        print("Score files saved.")

    return ngr_list, initial_scores

def normalize_matrix(W, file_name_prefix, normalization=''):
    # pickle file name
    matrix_pickle_file_name = file_name_prefix+'_'+str(W.shape[0])+'_genes_normalized'
    
    if len(normalization) == 0:
        print('\nMatrix normalization method not provided; default is personalized pageRank: W = D^-1/2 A D^-1/2')
    else:
        print('Matrix normalization method: {0}'.format(normalization))
        matrix_pickle_file_name = matrix_pickle_file_name + '_'+normalization
        
    matrix_pickle_file_name = matrix_pickle_file_name + '.pkl'

    if path.exists(matrix_pickle_file_name):
        W = pickle.load(open(matrix_pickle_file_name, 'rb'))
        print('\nNormalized matrix reading done. Existing pickle loaded.')
    else: # normalize matrix and dump the pickle
        D = np.diag(np.sum(W, axis=0)) # D is diagonal degree matrix

        if normalization == 'diffusion_kernel':
            W = D - W
        elif normalization in ('insulated_diffusion', 'electric_network'):
            W = np.matmul(W, np.linalg.inv(D))
        else: # personalized_pagerank (default)
            D12 = np.linalg.inv(np.sqrt(D))
            W = np.matmul(np.matmul(D12, W), D12)

        column_sums = np.sum(W, axis=0)
        if sum(column_sums > (1+0.05)) > 0: # sanity check (and application) for column-stochasticity of normalized matrix
            print('Matrix is not column-stochastic. One or more columns have sums > 1+epsilon.')
            #W = W / column_sums
            #print('Column-normalization performed.')
        
        with open(matrix_pickle_file_name, 'wb') as f:
            pickle.dump(W, f)

        print('\nMatrix normalization done. Pickle dumped.')

    print(datetime.datetime.now())
    return W


## Helper functions to evaluation
# sorts a structured numpy array (called list in this script) based on values in column colname
def sort_structured_array(input_array, colname='score', decreasing=True):
    sorted_order = np.argsort(input_array[colname])
    
    if decreasing:
        sorted_order = np.flip(sorted_order)
        
    return input_array[sorted_order]

# sets second ranked list as reference (with its ordered items ascending as 0...N, then
# maps the first ranked list based on these 0..N indices of the second; used in tau, manhattan distance, etc.
def rank_ranked_lists_relatively(input_list1, input_list2):
    map1 = {}
    for i in range(input_list1.shape[0]):
        map1[input_list1[i]] = i
        
    # loop over the second list and build vector to calculate euclidean space
    ranks1 = np.zeros((input_list1.shape[0], ))
    ranks2 = np.zeros((input_list2.shape[0], ))
    for i in range(input_list1.shape[0]):
        ranks2[i] = i
        ranks1[i] = map1[input_list2[i]] # rank of gene being iterated over in second list, in the first (gold standard) list as saved in the dictionary

    return ranks1, ranks2
