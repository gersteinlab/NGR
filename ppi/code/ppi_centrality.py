import sys, os, argparse
import networkx as nx
import pickle

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-k", help="number of sample nodes to calculate betweenness centrality")
    parser.add_argument("-nf", help="network filename")
    parser.add_argument("-nn", help="network name")
    
    args = parser.parse_args()

    if args.k:
        global k_nodes
        k_nodes = int(args.k)
        print("Parameter k set to {0}".format(k_nodes))
    
    if args.nf:
        global network_filename
        network_filename = args.nf
        print("Network filename set to {0}".format(network_filename))
        
    if args.nn:
        global network_name
        network_name = args.nn
        print("Network name set to {0}".format(network_name))
        
k_nodes = None # default is None and hence all nodes are used
network_filename = "../HumanNetv2_FN_converted.txt"
network_name = "HumanNetv2"

def get_filename(network_name, k_nodes, type='pickle', extension='pkl'):
    if k_nodes is None:
        filename_body = network_name
    else:
        filename_body = network_name+"_"+str(k_nodes)
    
    if type == 'pickle':
        filename = "pickle_bc_dict_"+filename_body
    else:
        filename = filename_body+"_"+type+"_values"
    
    filename = filename + "." + extension
    return filename

# returns a dictionary with gene names and their centrality values
def calculate_ppi_centrality_values(network_name, network_filename="", k_nodes=None):
    pickle_filename = get_filename(network_name, k_nodes, type='pickle', extension='pkl')
    if os.path.exists(pickle_filename):
        f = open(pickle_filename, "rb")
        bc = pickle.load(f)
        return bc    
    else: 
        G = nx.read_weighted_edgelist(network_filename)
        bc = nx.betweenness_centrality(G, k=k_nodes, normalized=True, weight='weight', endpoints=False, seed=18)
        print(bc)
        
        # pickle dictionary
        f = open(pickle_filename, "wb")
        pickle.dump(bc, f)
        print("Dictionaty pickle created.")
        
        return bc
    
def write_dict_to_CSV(input_dict, filename, header=""):

    dict_text=""
    if header != "":
        dict_text = header + "\n"
    
    for key, value in input_dict.items():
        dict_text = dict_text + (str(key)+","+str(value)) + "\n"
    
    try:
        output_file = open(filename, "w")
        output_file.write(dict_text)
        output_file.close()
        print("File {} written successfully.".format(filename))
    except:
        print("Error writing file {}.".format(filename))

if __name__ == "__main__":
    parse_arguments()
    bc = calculate_ppi_centrality_values(network_name, network_filename, k_nodes)
    
    results_filename = get_filename(network_name, k_nodes, type='bc', extension='csv')
    write_dict_to_CSV(bc, results_filename, header="gene,centrality")
    
