import numpy as np
import scipy.stats as st
from scipy.spatial import distance
from scipy.stats import weightedtau

import helper_functions as hp

# initial_scores is a structured array with two columns: 'gene' and 'score'
# result_genes and result_scores are the lists of genes and resulting scores to be evaluated
# gold_standard_file referes to csv file with two columns: 'gene' and 'score'
# comparison measure options: euclidean (Euclidean distance between individual genes w.r.t. gold standard)
def evaluate_results(ranked_lists, gold_standard_file, list_names, comp_measure='RBO', id_colname='gene', value_colname='score'):
    #  read gold standard structured array
    gold_standard_list = np.genfromtxt(gold_standard_file, dtype='U25,f8', names=('gene', 'score'), skip_header=1, delimiter=',')
    
    # evaluate each of the ranked lists against the gold standard
    for l in range(len(list_names)):        

        common_ids, gold_standard_inds, ranked_list_inds = np.intersect1d(gold_standard_list[id_colname], ranked_lists[l][id_colname], return_indices=True)
        
        # sort lists w.r.t. score_colname
        ranked_list = hp.sort_structured_array(ranked_lists[l][ranked_list_inds], decreasing=True)
        gold_standard_list = hp.sort_structured_array(gold_standard_list[gold_standard_inds], decreasing=True)
        
        print("["+list_names[l]+"]:", end="")
        if 'dist' in comp_measure:
            distance_measure = 'euclidean_dist'
            if 'manhattan' in comp_measure:
                distance_measure = 'manhattan_dist' 

            score = distance_measure(ranked_list, gold_standard_list, distance_measure=distance_measure, id_colname=id_colname)
        elif comp_measure == 'RBO':
            score = RBO(ranked_list, gold_standard_list, id_colname=id_colname)
        elif comp_measure == 'AO':
            score = AO(ranked_list, gold_standard_list, depth=ranked_list.shape[0], id_colname=id_colname)
        elif comp_measure == 'weighted_tau':
            score, _ = weighted_tau(ranked_list, gold_standard_list, id_colname=id_colname)
            
        print(score)

# EVALUATION FUNCTIONS:
# All evaluation functions below assume lists are sorted in decreasing order; each input list is a structured
# array with two columns titled 'gene' and 'score' unless noted otherwise

# calculates Eucliden or Manhattan distance between ranks of elements in two ordered lists
def distance_measure(input_list1, input_list2, distance_measure='manhattan_dist', id_colname='gene'):
    map1 = {}
    for i in range(input_list1.shape[0]):
        map1[input_list1[id_colname][i]] = i
        
    # loop over the second list and build vector to calculate euclidean space
    ranks1 = np.zeros((input_list1.shape[0], ))
    ranks2 = np.zeros((input_list2.shape[0], ))
    for i in range(input_list1.shape[0]):
        ranks2[i] = i
        ranks1[i] = map1[input_list2[id_colname][i]] # rank of gene being iterated over in second list, in the first list as saved in the dictionary

    if 'manhattan' in distance_measure:
        distance_value = np.sum(abs(ranks1 - ranks2))
    else:
        distance_value = distance.euclidean(ranks1, ranks2)
        
    return distance_value

# Average overlap as per implementation by Ritesh Agrawal at https://github.com/ragrawal/measures/blob/master/measures/rankedlist/
# See: https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
def AO(input_list1, input_list2, depth = 10, id_colname='gene'):
    l1 = list(input_list1[id_colname])
    l2 = list(input_list2[id_colname])

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
def RBO(input_list1, input_list2, p = 0.98, id_colname='gene'):
    l1 = list(input_list1[id_colname])
    l2 = list(input_list2[id_colname])
    
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

# A wrapper function for weighted tau. 
# See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.weightedtau.html
def weighted_tau(input_list1, input_list2, id_colname='gene'):
    map1 = {}
    for i in range(input_list1.shape[0]):
        map1[input_list1[id_colname][i]] = i
        
    # loop over the second list and build vector to calculate euclidean space
    ranks1 = np.zeros((input_list1.shape[0], ))
    ranks2 = np.zeros((input_list2.shape[0], ))
    for i in range(input_list1.shape[0]):
        ranks2[i] = i
        ranks1[i] = map1[input_list2[id_colname][i]] # rank of gene being iterated over in second list, in the first (gold standard) list as saved in the dictionary

    return st.weightedtau(ranks1, ranks2, rank=False)
    
