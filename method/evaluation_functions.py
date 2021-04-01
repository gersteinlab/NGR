import numpy as np
import scipy.stats as st
from scipy.spatial import distance
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score, precision_score, recall_score, precision_recall_curve, auc
from scipy.stats import hypergeom, rankdata

import helper_functions as hp
from scipy.constants.codata import precision

import sys
import pickle

# initial_scores is a structured array with two columns: 'gene' and 'score'
# result_genes and result_scores are 2-D structured arrays: lists of genes and resulting scores to be evaluated
# gold_standard_file referes to csv file with two columns: 'gene' and 'score'
# N is the number of genes to be compared across lists
# comparison measure options: euclidean (Euclidean distance between individual genes w.r.t. gold standard)
def compare_lists(ranked_score_lists, gold_standard_file, list_names, cancer_type, variant_type, comp_measures=('AO'), id_colname='gene', value_colname='score'):
    header = ''.ljust(20)
    for cm in comp_measures:
        header += '{0: <20}'.format(cm)
    print(header, end='')

    # evaluate each of the ranked lists against the gold standard
    for l in range(len(list_names)):
        gold_standard_scores = np.genfromtxt(gold_standard_file, dtype='U25,f8', names=('gene', 'score'), skip_header=1, delimiter=',')

        # retain expressed genes only in gold standrd list
        if cancer_type != '':
            cancer_type_expressed_genes_file = '../tcga/code/results/expressed_genes/'+cancer_type.lower()+'_expressed_genes.csv'
            cancer_type_expressed_genes = np.loadtxt(cancer_type_expressed_genes_file, dtype='U25', skiprows=1)
            
            _, cancer_type_expressed_genes_inds, gold_standard_expressed_gene_inds = np.intersect1d(cancer_type_expressed_genes, gold_standard_scores[id_colname], return_indices=True)               
            cancer_type_expressed_genes = cancer_type_expressed_genes[cancer_type_expressed_genes_inds]       
            gold_standard_scores = gold_standard_scores[gold_standard_expressed_gene_inds]

            # In case AUC of whole PPI is needed: comment second line re: cancer_type_expressed_genes = ...m, and uncomment 
            # the block below to locate unexpressed PPI genes to the bottom of the gold standard list with 0 values

            #unexpressed_ranked_genes = np.setdiff1d(ranked_score_lists[l][id_colname], cancer_type_expressed_genes)
            #unexpressed_ranked_gene_scores = np.zeros(len(unexpressed_ranked_genes))
            #unexpressed_ranked_genes_arr = np.array(list(zip(unexpressed_ranked_genes, unexpressed_ranked_gene_scores)), dtype=gold_standard_scores.dtype)
            #gold_standard_scores = np.append(gold_standard_scores, unexpressed_ranked_genes_arr)
            #gold_standard_scores = np.sort(gold_standard_scores, order='gene')

        common_ids, gold_standard_inds, ranked_score_inds = np.intersect1d(gold_standard_scores[id_colname], ranked_score_lists[l][id_colname], return_indices=True)
        N = len(common_ids) # number of post-processed genes
        
        # Re-sort lists in decreasing order of score
        gold_standard_scores = np.flip(np.sort(gold_standard_scores[gold_standard_inds], order='score'))
        ranked_scores = np.flip(np.sort(ranked_score_lists[l][ranked_score_inds], order='score'))

        print('\n{0: <20}'.format(list_names[l]), end="")
        for cm in comp_measures:
            if 'dist' in cm:
                dist_measure = 'euclidean_dist'
                if 'manhattan' in cm:
                    dist_measure = 'manhattan_dist' 
                
                score = distance_measure(ranked_scores[id_colname], gold_standard_scores[id_colname], measure=dist_measure)
                score /= (ranked_scores.shape[0]*ranked_scores.shape[0]) # normalize distance
            elif cm == 'RBO':
                score = RBO(ranked_scores[id_colname], gold_standard_scores[id_colname])
            elif cm == 'AO':
                score = AO(ranked_scores[id_colname], gold_standard_scores[id_colname], depth=ranked_scores.shape[0])
            elif 'tau' in cm:
                if 'weighted' in cm:
                    score, _ = weighted_tau(ranked_scores[id_colname], gold_standard_scores[id_colname])
                else:
                    score, _ = tau(ranked_scores[id_colname], gold_standard_scores[id_colname])
            elif cm in ('F1', 'Accuracy', 'AUC'):
                rank_threshold = int(0.1 * N) # label 1 for ranks above threshold; 0 otherwise
                gold_standard_labels = np.concatenate((np.ones(rank_threshold), np.zeros(N-rank_threshold)))
                
                # ngr labels reordered so labels match and AUC is measured
                ngr_labels = np.concatenate((np.ones(rank_threshold), np.zeros(N-rank_threshold)))

                gold_standard_order, _ = hp.rank_ranked_lists_relatively(ranked_scores['gene'], gold_standard_scores['gene'])
                ngr_labels = ngr_labels[gold_standard_order]

                if cm == 'F1':
                    score = f1_score(gold_standard_labels, ngr_labels)
                elif cm == 'Accuracy':
                    score = accuracy_score(gold_standard_labels, ngr_labels)
                elif cm == 'Precision':
                    score = precision_score(gold_standard_labels, ngr_labels)
                elif cm == 'AUC':                
                    score = roc_auc_score(gold_standard_labels, ngr_labels)
                elif cm == 'AUPRC':
                    precision, recall, thresholds = precision_recall_curve(gold_standard_labels, ngr_labels)
                    score = auc(recall, precision)
            elif cm == 'Hypergeom_test': # hypergeometric enrichment test
                known_genes_file = '../gene_lists/known_genes/COSMIC_v90_membership.csv'
                known_genes = np.genfromtxt(known_genes_file, delimiter=',', skip_header=1, dtype='U25,i4', names=('gene','score'))

                S = int(0.1 * N) #int(0.025 * N) # sample size
 
                sample_common_known_genes = np.intersect1d(known_genes['gene'], ranked_scores[id_colname][0:S])
                x = len(sample_common_known_genes) # number of "successes" in the sample of size N

                population_common_known_genes = np.intersect1d(known_genes['gene'], ranked_scores[id_colname])
                n = len(population_common_known_genes) # number of "successes" in entire population

                pval = hypergeom.sf(x-1, N, n, S)
                score = 0
                if pval < 0.05:
                    score = 1

                print('{0}, {1}, {2}, {3}. P-val: {4}'.format(N, n, S, x, pval), end=" ")
                
            score = round(score, 3)
            print('{0: <20}'.format(score), end="")

## EVALUATION FUNCTIONS:

# MOBILITY LIST FUNCTIONS
# generated mobility lists across variant types, cancer types, ppis and ppi_levels
def generate_batch_gene_mobility_lists(uids_filename):
    matrices_dir = 'results/score_matrices/'
    ppi_matrix_dir = '../ppi/code/results/'
    genomics_matrix_dir = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'
    output_dir = 'results/mobility_lists/whole_lists/'

    uids_file = open(uids_filename, 'r')
    lines = uids_file.readlines()

    for l in lines:
        values = l.split(' ')
        variant_type = values[0]
        cancer_type = values[1]
        ppi_network = values[2]
        ppi_level = values[3]
        uid = values[4].strip()
        
        print('\n[{0}, {1}, {2}, {3} started. UID: {4}.]'.format(variant_type, cancer_type, ppi_network, ppi_level, uid))
                          
        output_pickle = matrices_dir+variant_type+'_'+cancer_type+'_'+ppi_network+'_lcc_'+cancer_type.lower()+'_'+uid+'_score_matrix.pkl'
        ngr_score_matrix = pickle.load(open(output_pickle, "rb"))

        ppi_matrix_prefix = ppi_matrix_dir+ppi_network+'_converted_matrix'
        genomics_matrix_prefix = genomics_matrix_dir+variant_type+'_'+cancer_type+'_matrix'
        
        ppi_matrix, ppi_matrix_index, extended_ppi_matrix_prefix, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = hp.get_ngr_inputs(genomics_matrix_prefix, ppi_matrix_prefix, genomics_matrix_transformed=True, largest_cc=True, ppi_cancer_type=cancer_type, matrix_normalization_method='insulated_diffusion')
        
        mobility_list = generate_gene_mobility_list(genomics_matrix, ngr_score_matrix, genomics_matrix_row_index)
        mobility_list_filename = output_dir+variant_type+'_'+cancer_type+'_'+ppi_network+'_lcc_'+cancer_type.lower()+'_mobility_list.csv'
        
        np.savetxt(mobility_list_filename, np.stack((mobility_list['gene'], mobility_list['mobility_score'], mobility_list['initial_rank'], mobility_list['ngr_rank']), axis=1), fmt='%s', delimiter=',')

        print('\n[{0}, {1}, {2}, {3} mobility list saved in {4}]'.format(variant_type, cancer_type, ppi_network, ppi_level, mobility_list_filename))

# generates lists sorted in increasing order of mobility: difference between initial and post-diffusion rank
# values are -ve (for upward mobility) and +ve (downward); upward 'passenger' genes are of more interest
def generate_gene_mobility_list(initial_matrix, ngr_matrix, gene_index):
    initial_ranks = rankdata(-np.mean(initial_matrix, axis=1), method='min')
    ngr_ranks = rankdata(-np.mean(ngr_matrix, axis=1), method='min')
    
    mobility_score = initial_ranks - ngr_ranks # mobility is the different in ranks between initial and ngr list
    
    mobility_list = np.array(list(zip(gene_index, mobility_score, initial_ranks, ngr_ranks)), dtype=[('gene', 'U25'), ('mobility_score', 'i4'), ('initial_rank', 'i4'), ('ngr_rank', 'i4')])
    mobility_list = np.flip(np.sort(mobility_list, order=['mobility_score']))
    
    return mobility_list

# LIST COMPARISON FUNCTIONS
# All evaluation functions below assume lists are sorted in decreasing order; each input list is 1-D (i.e. a column from the structured array used in compare_lists() above
# In these functions, a list is 1-D

def generate_batch_hypergemetric_pvalues(uids_filename):
    matrices_dir = 'results/score_matrices/'
    ppi_matrix_dir = '../ppi/code/results/'
    genomics_matrix_dir = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'
    gold_standard_file = '../gene_lists/code/results/gold_standard_lists/cancerMine_both_depmap_both_values_mean_scores_breast_tissue.csv' # used only to run the function successfully

    uids_file = open(uids_filename, 'r')
    lines = uids_file.readlines()

    for l in lines:
        values = l.split(' ')
        variant_type = values[0]
        cancer_type = values[1]
        ppi_network = values[2]
        ppi_level = values[3]
        uid = values[4].strip()
        
        print('\nWorking on {0}, {1}, {2}, {3}: {4}'.format(variant_type, cancer_type, ppi_network, ppi_level, uid))
             
        output_pickle = matrices_dir+variant_type+'_'+cancer_type+'_'+ppi_network+'_lcc_'+cancer_type.lower()+'_'+uid+'_score_matrix.pkl'
        ngr_score_matrix = pickle.load(open(output_pickle, "rb"))

        ppi_matrix_prefix = ppi_matrix_dir+ppi_network+'_converted_matrix'
        genomics_matrix_prefix = genomics_matrix_dir+variant_type+'_'+cancer_type+'_matrix'
        
        ppi_matrix, ppi_matrix_index, extended_ppi_matrix_prefix, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = hp.get_ngr_inputs(genomics_matrix_prefix, ppi_matrix_prefix, genomics_matrix_transformed=True, largest_cc=True, ppi_cancer_type=cancer_type, matrix_normalization_method='insulated_diffusion')

        # ngr scores in a structured array
        ngr_scores = hp.process_ngr_results(ngr_score_matrix, ppi_matrix_index, save_files=False, uid=uid) 

        compare_lists((ngr_scores, ngr_scores), gold_standard_file, ('Tissue-PPI', 'Tissue-PPI-Redundant'), cancer_type=cancer_type, variant_type=variant_type, comp_measures=('Accuracy', 'Hypergeom_test'))
        
        
# calculates Eucliden or Manhattan distance between ranks of elements in two ordered lists
def distance_measure(input_list1, input_list2, measure='manhattan_dist'):
    ranks1, ranks2 = hp.rank_ranked_lists_relatively(input_list1, input_list2)

    if 'manhattan' in measure:
        distance_value = np.sum(abs(ranks1 - ranks2))
    else:
        distance_value = distance.euclidean(ranks1, ranks2)

    return distance_value

# Average overlap as per implementation by Ritesh Agrawal at https://github.com/ragrawal/measures/blob/master/measures/rankedlist/
# See: https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
def AO(input_list1, input_list2, depth = 10):
    l1 = list(input_list1)
    l2 = list(input_list2)

    if l1 == None: l1 = []
    if l2 == None: l2 = []

    sl, ll = sorted([(len(l1), l1),(len(l2),l2)])
    s, S = sl  # s = length of smaller list, S = Smaller List
    l, L = ll  # l = length of longer list, L = Longer list
    #sanity check
    if s == 0: return 0
    depth = depth if depth < l else l
    
    # Calculate fraction of overlap from rank  at ranks 1 through depth
    # (the longer of the two lists)
    ss = set([])
    ls = set([])
    overlap = {0: 0}  # overlap holds number of common elements at depth d 
    sum1 = 0.0  

    for i in range(depth):
        # get elements from the two list
        x = L[i]
        y = S[i] if i < s else None
        depth = i+1
        # if the two elements are same, then we don't need
        # to them to the list and just increment the 
        if x == y: 
            overlap[depth] = overlap[i] + 2
        #else add items to the two list
        else:
            ls.add(x)
            if y != None: ss.add(y)
            overlap[depth] = overlap[i] + (2 if x in ss else 0) + (2 if y in ls else 0) 
        sum1 = sum1 + float(overlap[depth])/(len(S[0:depth]) + depth)

    ao_score = sum1/depth
    return ao_score

# Rank Bias Overlap score as per implementation by Ritesh Agrawal at https://github.com/ragrawal/measures/blob/master/measures/rankedlist/
# See: https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
def RBO(input_list1, input_list2, p = 0.98):
    l1 = list(input_list1)
    l2 = list(input_list2)
    
    if l1 == None: l1 = []
    if l2 == None: l2 = []
    
    sl,ll = sorted([(len(l1), l1),(len(l2),l2)])
    s, S = sl
    l, L = ll
    if s == 0: return 0

    # Calculate the overlaps at ranks 1 through l 
    # (the longer of the two lists)
    ss = set([]) # contains elements from the smaller list till depth i
    ls = set([]) # contains elements from the longer list till depth i
    x_d = {0: 0}
    sum1 = 0.0
    for i in range(l):
        x = L[i]
        y = S[i] if i < s else None
        d = i + 1
        
        # if two elements are same then 
        # we don't need to add to either of the set
        if x == y: 
            x_d[d] = x_d[d-1] + 1.0
        # else add items to respective list
        # and calculate overlap
        else: 
            ls.add(x) 
            if y != None: ss.add(y)
            x_d[d] = x_d[d-1] + (1.0 if x in ss else 0.0) + (1.0 if y in ls else 0.0)     
        #calculate average overlap
        sum1 += x_d[d]/d * pow(p, d)
        
    sum2 = 0.0
    for i in range(l-s):
        d = s+i+1
        sum2 += x_d[d]*(d-s)/(d*s)*pow(p,d)

    sum3 = ((x_d[l]-x_d[s])/l+x_d[s]/s)*pow(p,l)

    # Equation 32
    rbo_score = (1-p)/p*(sum1+sum2)+sum3
    return rbo_score

# A wrapper function for weighted Kendall tau measure 
# See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weightedtau.html
def weighted_tau(input_list1, input_list2):
    ranks1, ranks2 = hp.rank_ranked_lists_relatively(input_list1, input_list2)
    return st.weightedtau(ranks2, ranks1, rank=False)

# A wrapper function for Kendall tau measure
# See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html
def tau(input_list1, input_list2):
    ranks1, ranks2 = hp.rank_ranked_lists_relatively(input_list1, input_list2)
    return st.kendalltau(ranks2, ranks1)
