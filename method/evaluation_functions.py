import helper_functions as hp
import numpy as np

# initial_scores is a structured array with two columns: 'gene' and 'score'
# result_genes and result_scores are the lists of genes and resulting scores to be evaluated
# gold_standard_file referes to csv file with two columns: 'gene' and 'score'
# comparison measure options: euclidean (Euclidean distance between individual genes w.r.t. gold standard)
def evaluate_results(initial_scores, ngr_list, gold_standard_file, comp_measure='euclidean', save_files=True):
    #  read gold standard structured array
    gold_standard_list = np.genfromtxt(gold_standard_file, dtype='U25,f8', names=('gene', 'score'), skip_header=1, delimiter=',')

    # compare lists: can add more lists to the input tuple below; currently only two lists are used as needed
    if comp_measure in ['euclidean']:
        compare_ranked_lists((initial_scores, ngr_list), ('initial', 'ngr'), 'gene', 'score', gold_standard_list, comp_measure=comp_measure)
    # else: to allow for scalability, might for example later decide to evalaute based on top k% (known) genes in resulting list rather than order comparison with initial scores

# ranked lists is a tuple of structured arrays; same set of elements across lists is assumed (e.g. same number and IDs of genes)
# names are descriptors of each ndarray lists used to report final results
# value_colname is the name of column based on which sorting/ranking is done 
# id_colname is the name of column with element ids based on which comparison is done (e.g. gene names)
def compare_ranked_lists(ranked_lists, names, id_colname, value_colname, gold_standard_list, comp_measure='euclidean'):

    # find intersecting elements
    for l in range(len(names)):
        print("\n\n"+names[l])
        
        common_ids, gold_standard_inds, ranked_list_inds = np.intersect1d(gold_standard_list[id_colname], ranked_lists[l][id_colname], return_indices=True)

        # sort gold standard w.r.t. score_colname
        gold_standard_list = hp.sort_structured_array(gold_standard_list[gold_standard_inds], value_colname, decreasing=True)
        ranked_list = hp.sort_structured_array(ranked_lists[l][ranked_list_inds], value_colname, decreasing=True)

        hp.calculate_euclidean_distance_basic(ranked_list, gold_standard_list, id_colname)

    