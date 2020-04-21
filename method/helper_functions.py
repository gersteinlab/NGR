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

# reads input matrix and its (row and col) indices
# In matrix file: rows correpond to samples, cols to genes; We transpose in NGR's numpy frame (note the .T below),
# where rows end up representing (N) genes and cols (P) sample vectors. 
def read_matrix(file_name_prefix, sorted=True, target='both', selected_entities={'rows':[], 'cols':[]}):
    matrix_pickle_file_name = file_name_prefix+'.pkl'
    matrix_pickle_row_index_file_name = file_name_prefix+'_gene_index.pkl'
    matrix_pickle_col_index_file_name = file_name_prefix+'_barcode_index.pkl'

    if path.exists(matrix_pickle_file_name) and path.exists(matrix_pickle_row_index_file_name) and path.exists(matrix_pickle_col_index_file_name):
        matrix = pickle.load(open(matrix_pickle_file_name, 'rb'))
        matrix_row_index = pickle.load(open(matrix_pickle_row_index_file_name, 'rb'))
        matrix_col_index = pickle.load(open(matrix_pickle_col_index_file_name, 'rb'))
        print('\nMatrix and indices reading done. Existing pickles loaded.')
    else:
        matrix_file_name = file_name_prefix+'.txt'
        matrix_row_index_file_name = file_name_prefix+'_gene_index.csv'
        matrix_col_index_file_name = file_name_prefix+'_barcode_index.csv'
        
        matrix = np.loadtxt(matrix_file_name, dtype='f').T # pickle matrix
        matrix_row_index = np.loadtxt(matrix_row_index_file_name, dtype='U25') # pickle matrix index
        matrix_col_index = np.loadtxt(matrix_col_index_file_name, dtype='U25') # pickle matrix index
        
        if sorted:
            matrix, matrix_row_index, matrix_col_index = sort_matrix(matrix, matrix_row_index, matrix_col_index, target=target)
        
        if (len(selected_entities['rows']) > 0 or len(selected_entities['cols']) > 0):
            if sorted:
                matrix, matrix_row_index, matrix_col_index = select_from_matrix(matrix, matrix_row_index, matrix_col_index, selected_entities, target=target)
            else:
                raise Exception('Warning: Target matrix index(ces) not sorted. Cannot select entities.'.format(target))

        pickle.dump(matrix, open(matrix_pickle_file_name, 'wb'))
        pickle.dump(matrix_row_index, open(matrix_pickle_row_index_file_name, 'wb'))
        pickle.dump(matrix_col_index, open(matrix_pickle_col_index_file_name, 'wb'))    
        print('\nMatrix and indices reading done. Pickles saved.')
    
    print('{0}\n'.format(datetime.datetime.now()))
    return matrix, matrix_row_index, matrix_col_index

# sort matrix based on row or col names (or both)
def sort_matrix(input_matrix, input_matrix_row_index, input_matrix_col_index, target='both'):
    matrix = input_matrix
    matrix_row_index = input_matrix_row_index
    matrix_col_index = input_matrix_col_index

    if target in ('rows', 'both'): # sort by row name
        sorted_order = np.argsort(matrix_row_index)
        matrix_row_index = matrix_row_index[sorted_order]
        matrix = matrix[sorted_order, :]
    
    if target in ('cols', 'both'):
        sorted_order = np.argsort(matrix_col_index)
        matrix_col_index = matrix_col_index[sorted_order]
        matrix = matrix[:, sorted_order]

    print('Matrix sorting of {0} done.'.format(target))

    return matrix, matrix_row_index, matrix_col_index

def select_from_matrix(input_matrix, input_matrix_row_index, input_matrix_col_index, selected_entities, target='both'):
    matrix = input_matrix
    matrix_row_index = input_matrix_row_index
    matrix_col_index = input_matrix_col_index

    if target in ('rows', 'both') and len(selected_entities['rows']) > 0:
        selected_entities['rows'] = np.sort(selected_entities['rows'])                                  
        common_entities = np.intersect1d(matrix_row_index, selected_entities['rows'])
        common_entities_in_matrix_index = np.searchsorted(matrix_row_index, common_entities) # indices of common entities in the matrix
        matrix = matrix[common_entities_in_matrix_index, :] # select from matrix rows 
        matrix_row_index = matrix_row_index[common_entities_in_matrix_index] # select from index
    if target in ('cols', 'both') and len(selected_entities['cols']) > 0:
        selected_entities['cols'] = np.sort(selected_entities['cols'])
        common_entities = np.intersect1d(matrix_col_index, selected_entities['cols'])
        common_entities_in_matrix_index = np.searchsorted(matrix_col_index, common_entities) # indices of common entities in the matrix
        matrix = matrix[:, common_entities_in_matrix_index]
        matrix_col_index = matrix_col_index[common_entities_in_matrix_index]

    print('\nMatrix and {0} index selection done.'.format(target))
    print('Final matrix dimensions: ' + str(matrix.shape))
    print('{0}\n'.format(datetime.datetime.now()))
    
    return matrix, matrix_row_index, matrix_col_index

# reads NxN ppi matrix and its index
def read_square_matrix(file_name_prefix, sorted=True, largest_cc=True, selected_entities=[]):
    matrix_pickle_file_name = file_name_prefix+'.pkl'
    matrix_pickle_index_file_name = file_name_prefix+'_index.pkl'
    if largest_cc:
        matrix_pickle_file_name = file_name_prefix+'_lcc.pkl'
        matrix_pickle_index_file_name = file_name_prefix+'_lcc_index.pkl'

    if path.exists(matrix_pickle_file_name) and path.exists(matrix_pickle_index_file_name):
        matrix = pickle.load(open(matrix_pickle_file_name, 'rb'))
        matrix_index = pickle.load(open(matrix_pickle_index_file_name, 'rb'))
        print('\nSquare matrix and index reading done. Existing pickles loaded.')
    else:
        matrix_file_name = file_name_prefix+'.txt'
        matrix_index_file_name = file_name_prefix+'_index.txt'
        
        matrix = np.loadtxt(matrix_file_name, dtype='f') # pickle matrix
        matrix_index = np.loadtxt(matrix_index_file_name, dtype='U25') # pickle matrix index
        
        if largest_cc:
            matrix, matrix_index = get_largest_connected_component(matrix, matrix_index)
        
        if sorted:
            matrix, matrix_index = sort_square_matrix(matrix, matrix_index)
        
        if (len(selected_entities) > 0):
            if sorted:
                matrix, matrix_index = select_from_square_matrix(matrix, matrix_index, selected_entities)
            else:
                raise Exception('Warning: Square matrix index not sorted. Cannot select entities.'.format(target))

        pickle.dump(matrix, open(matrix_pickle_file_name, 'wb'))
        pickle.dump(matrix_index, open(matrix_pickle_index_file_name, 'wb'))
        print('\nSquare matrix and index reading done. Pickles saved.')
    
    print('{0}\n'.format(datetime.datetime.now()))
    return matrix, matrix_index
    
# sorts NxN matrix according to given index that represents both cols and rows
def sort_square_matrix(input_matrix, input_matrix_index):
    matrix = input_matrix
    matrix_index = input_matrix_index
    
    sorted_order = np.argsort(matrix_index)
    matrix = matrix[sorted_order, :]
    matrix = matrix[:, sorted_order]
    matrix_index = matrix_index[sorted_order]
    print('Matrix row-col sorting done.')
 
    return matrix, matrix_index

# selects from NxN matrix selected_entities that represent both rows and cols
def select_from_square_matrix(input_matrix, input_matrix_index, selected_entities):
    matrix = input_matrix
    matrix_index = input_matrix_index

    selected_indices = np.sort(selected_entities)
    common_entities = np.intersect1d(matrix_index, selected_entities)
    common_entities_in_matrix_index = np.searchsorted(matrix_index, common_entities) # indices of common entities in the matrix
    print(common_entities_in_matrix_index)
    
    matrix = matrix[:, common_entities_in_matrix_index] # reorder matrix rows and cols
    matrix = matrix[common_entities_in_matrix_index, :]

    matrix_index = matrix_index[common_entities_in_matrix_index] # reorder index

    print('\nSquare matrix and index selection done.')
    print('Final square matrix dimensions: ' + str(matrix.shape))
    print('{0}\n'.format(datetime.datetime.now()))
    
    return matrix, matrix_index

# returns an updated matrix including new edges from input edge list
# IMPORTANT NOTE: new edges are by default directed; ppi's are originally undirected (but are rendered directed by definition once new directed edges are added)
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
                        
    print('New edges have been added from {0}'.format(input_edgelist_file))
    
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
    
def get_ngr_inputs(genomics_matrix_file_prefix, ppi_matrix_file_prefix, sorted=True, largest_cc=True, matrix_normalization_method='insulated_diffusion'):
    ppi_matrix, ppi_matrix_index = read_square_matrix(ppi_matrix_file_prefix, sorted=sorted, largest_cc=largest_cc) # NxN
    
    genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = read_matrix(genomics_matrix_file_prefix, sorted=sorted, target='both') # N'xP, N' << N
    genomics_matrix, genomics_matrix_row_index = expand_genomics_matrix(ppi_matrix, ppi_matrix_index, genomics_matrix, genomics_matrix_row_index) # NxP
    
    if matrix_normalization_method != 'none': # ppi normalization
        if largest_cc:
            ppi_matrix = normalize_matrix(ppi_matrix, ppi_matrix_file_prefix+'_lcc', normalization=matrix_normalization_method)
        else:
            ppi_matrix = normalize_matrix(ppi_matrix, ppi_matrix_file_prefix, normalization=matrix_normalization_method)

    return ppi_matrix, ppi_matrix_index, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index

# expands the input genomics matrix to include all (and only) genes in the ppi; genes in genomics matrix
# not in ppi are removed; genes in ppi not in matrix added with zero rows across patients to allow for diffusion
# genomics matrix row index updated accordingly
def expand_genomics_matrix(ppi_matrix, ppi_matrix_index, genomics_matrix, genomics_matrix_row_index):
    common_entities = np.intersect1d(ppi_matrix_index, genomics_matrix_row_index)
    common_entities_indices_in_genomics_matrix_row_index = np.searchsorted(genomics_matrix_row_index, common_entities)
    common_entities_indices_in_ppi_matrix_index = np.searchsorted(ppi_matrix_index, common_entities)
    
    expanded_genomics_matrix = np.zeros((ppi_matrix.shape[0], genomics_matrix.shape[1]))
    expanded_genomics_matrix[common_entities_indices_in_ppi_matrix_index, :] = genomics_matrix[common_entities_indices_in_genomics_matrix_row_index, :]

    expanded_genomics_matrix_row_index = ppi_matrix_index
    
    print('Input matrix expansion done. Expanded matrix dimensions: {0}.'.format(expanded_genomics_matrix.shape))
    return expanded_genomics_matrix, expanded_genomics_matrix_row_index


# initial_scores is a structured array with two columns: 'gene' and 'score'
# ngr_scores is a numpy arrays of scores
# processes the output of NGR and returns two structured arrays to facilitate evaluation; each array has two columns, 'gene' and 'score,' and its rows are sorted in decreasing order of score
def process_ngr_results(ngr_score_matrix, ngr_genes, save_files=True, uid='', output_dir=''):
    # create structures array and sort result scores in decreasing order 
    dtype = [('gene', 'U25'), ('score', 'f')]
    ngr_scores = np.array(list(zip(ngr_genes, np.mean(ngr_score_matrix, axis=1))), dtype=dtype)
    ngr_scores = np.flip(np.sort(ngr_scores, order='score'))
    
    if save_files:
        if uid != '':
            uid = str(uid)+'_'
            
        np.savetxt(output_dir+uid+'ngr_scores.csv', np.stack((ngr_scores['gene'], ngr_scores['score']), axis=1), fmt='%s', delimiter=',')
        print("Score files saved.")

    return ngr_scores

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
        
        pickle.dump(W, open(matrix_pickle_file_name, 'wb'))
        print('\nMatrix normalization done. Pickle dumped.')

    print(datetime.datetime.now())
    return W

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

# sets second ranked list as reference (with its ordered entities ascending as 0...N, then
# maps the first ranked list based on these 0..N indices of the second; used in tau, manhattan distance, etc.
def rank_ranked_lists_relatively(input_list1, input_list2):
    map1 = {}
    for i in range(input_list1.shape[0]):
        map1[input_list1[i]] = i
        
    # loop over the second list and build vector to calculate euclidean space
    ranks1 = np.zeros((input_list1.shape[0], ), dtype=int)
    ranks2 = np.zeros((input_list2.shape[0], ), dtype=int)
    for i in range(input_list1.shape[0]):
        ranks2[i] = i
        ranks1[i] = map1[input_list2[i]] # rank of gene being iterated over in second list, in the first (gold standard) list as saved in the dictionary

    return ranks1, ranks2
