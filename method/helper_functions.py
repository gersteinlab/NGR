import numpy as np
from scipy.special import softmax

# Y is the vector of initial scores
# W is the weight matrix, i.e. PPI network
# W_index is the ordered list of genes as they appear in cols and rows
def ngr(Y, W, W_index):
    Y = Y.reshape(Y.shape[0], 1)
    S = np.ndarray.copy(Y)
    print(Y.shape)
    print(S.shape)
    
    iterations = 4
    i = 0
    while i < iterations:
        S = np.repeat(S, S.shape[0], axis=1)
        Sm = (W > 0) * S # Step 1: W-based ReLU on S 
        Sm = np.apply_along_axis(softmax, axis=1, arr=Sm) # Step 2: row-wise softmax
        Sm = np.transpose(Sm) # Step 3: Sm = Sm^T
        S = np.matmul(W, Sm) # Step 4: S = W.Sm
        S = np.diag(S).reshape(S.shape[0], 1) # Step 5: updated scores vector S = diagonal of mresulting matrix S

        i += 1
        print('Iteration ' +str(i)+ ' complete.')
        
    # start here: Helpful tips if needed: implement it in place and with 0-based softmax
    
    return S
    
# score files have two columns: gene names and scores
def read_scores(file_name, sorted=True):
    scores = np.loadtxt(file_name, delimiter=',', skiprows=1, 
                        dtype= {'names': ('gene', 'score'), 'formats': ('U20', 'f')})
    
    if sorted:
        sorted_order = np.argsort(scores['gene'])
        scores = scores[sorted_order]

    return scores

def read_matrix(file_name_prefix, sorted=True, selected_genes=[]):
    matrix = np.loadtxt(file_name_prefix+".txt", dtype='f')
    matrix_index = np.loadtxt(file_name_prefix+"_index.txt", dtype='U25')
    print("Matrix and index reading done.")
    
    if sorted:
        sorted_order = np.argsort(matrix_index)
        matrix = matrix[sorted_order, :]
        matrix = matrix[:, sorted_order]
        print("Matrix row-col sorting done.")
 
    if (len(selected_genes) > 0): # select genes in rows and columns
        #if(np.diff(selected_genes) >= 0 or np.diff(matrix_index) >= 0):
        #    raise Exception("Warning: selected_genes and/or matrix index are not sorted.")

        common_genes = np.intersect1d(matrix_index, selected_genes)
        common_genes_in_matrix_index = np.searchsorted(matrix_index, common_genes) # indices of common genes in the matrix

        print(selected_genes)
        print(common_genes_in_matrix_index)
        
        matrix = matrix[:, common_genes_in_matrix_index] # reorder matrix rows and cols
        matrix = matrix[common_genes_in_matrix_index, :]
        
        matrix_index = matrix_index[common_genes_in_matrix_index] # reorder index

    return matrix, matrix_index

# scaling divides all values in the matrix by the matrix-wide maximum to shrink values to [0, 1]
def normalize_matrix(W, scaling=True, normalization='insulated_diffusion'):
    
    if len(normalization) == 0:
        print('Normalization method not provided; default is personalized pageRank: W = D^-1/2 A D^-1/2')
    
    if scaling:
        W /= np.max(W)
        
    D = np.diag(np.sum((W > 0), axis=1)) # D is diagonal degree matrix
    if normalization == 'diffusion_kernel':
        W = D - W
    elif normalization in ('insulated_diffusion', 'electric_network'):
        W = np.matmul(W, np.linalg.inv(D))
    else:
        # personalized_pagerank (default)
        D12 = np.linalg.inv(np.sqrt(D))
        W = np.matmul(np.matmul(W, D12), D12)

    return W





