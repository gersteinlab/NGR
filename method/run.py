import argparse

import ngr
import helper_functions as hp
import evaluation_functions as eval

import numpy as np
from scipy.spatial import distance

import datetime
import uuid
import pickle

parser = argparse.ArgumentParser(description="Argument Parser")
parser.add_argument("--ct", help="cancer type")
parser.add_argument("--ppi", help="ppi network")
args = parser.parse_args()

cancer_type = 'LUSC' 
ppi_network='STRING'

if args.ct:
    cancer_type = args.ct
    print("Cancer type set to {}".format(cancer_type))

if args.ppi:
    ppi_network = args.ppi
    print("PPI network set to {}".format(ppi_network))

variant_type = 'somatic_MC3';  output_dir='results/score_matrices/'; output_filename_suffix = 'lcc_'+cancer_type.lower()    
ppi_matrix_prefix = '../ppi/code/results/'+ppi_network+'_converted_matrix'
genomics_matrix_prefix = '../tcga/cortex_data/variant_data/annovar/results/matrix_results/'+variant_type+'_'+cancer_type+'_matrix'

def main():
    uid = str(uuid.uuid4())[0:8] # unique id
    print('UID: ' +str(uid))

    # get data (ppi_matrix and scores)
    ppi_matrix, ppi_matrix_index, extended_ppi_matrix_prefix, genomics_matrix, genomics_matrix_row_index, genomics_matrix_col_index = hp.get_ngr_inputs(genomics_matrix_prefix, ppi_matrix_prefix, genomics_matrix_transformed=True, largest_cc=True, ppi_cancer_type=cancer_type, matrix_normalization_method='insulated_diffusion')
    print('PPI matrix: {0}, Input matrix: {1}.\n'.format(ppi_matrix.shape, genomics_matrix.shape))

    '''
    # (optional block) add new pathway edges and normalize ppi matrix
    #ppi_matrix[ppi_matrix > 0] = 1 # binarize matrix before adding edges
    #new_edges_file = '../ppi/curated_pathway_edges/SanchezVega-etal-pathway_edges.txt'
    #new_edges_suffix = 'sanchez_etal_pathway_edges'
    #output_filename_suffix += '_w_'+new_edges_suffix
     
    #ppi_matrix = hp.add_new_edges(ppi_matrix, ppi_matrix_index, new_edges_file)
    #ppi_matrix = hp.normalize_matrix(ppi_matrix, ppi_matrix_prefix+'_'+output_filename_suffix, normalization='insulated_diffusion') # normalize matrix after adding new edges
    '''
    
    # run method
    ngr_score_matrix = ngr.run(X=genomics_matrix, W=ppi_matrix)
    pickle.dump(ngr_score_matrix, open('_'.join([output_dir+variant_type, cancer_type, ppi_network, output_filename_suffix, uid, 'score_matrix.pkl']), 'wb'))
 
    '''
    ## generate removed gene indices lists (i.e. unexpressed ones in whole PPIs
    # once all results are generated; modify code block accordingly to include both cancer type and ppi name

    cancer_type_expressed_genes_file = '../tcga/code/results/expressed_genes/'+cancer_type.lower()+'_expressed_genes.csv'
    cancer_type_expressed_genes = np.loadtxt(cancer_type_expressed_genes_file, dtype='U25', skiprows=1)
    _, cancer_type_expressed_genes_inds, gold_standard_expressed_gene_inds = np.intersect1d(cancer_type_expressed_genes, gold_standard_scores[id_colname], return_indices=True)

    # In case finding ranks of removed (unexpressed) genes is needed to compare whole vs tissue specific PPI results, uncomment block below; figure to be created using R
    if(len(np.setdiff1d(ranked_score_lists[l][id_colname], cancer_type_expressed_genes)) > 0):
        removed_gene_indices = hp.find_removed_gene_indices(ranked_score_lists[l][id_colname], cancer_type_expressed_genes)
        removed_gene_indices_filename = 'results/removed_gene_indices/'+cancer_type.lower()+'_removed_gene_indices.csv'
        with open(removed_gene_indices_filename, 'w') as filehandle:
            for index in removed_gene_indices:
                filehandle.write('%s\n' % index)
                
    '''

main()
