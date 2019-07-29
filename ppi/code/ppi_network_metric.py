import sys, os, argparse
import networkx as nx
import pickle

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-k", help="number of sample nodes to calculate betweenness centrality")
    parser.add_argument("-nf", help="network filename")
    parser.add_argument("-nn", help="network name")
    parser.add_argument("-mt", help="metric")
    parser.add_argument("-wt", help="edge weight flag")
    
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
    
    if args.mt:
        global metric
        metric = args.mt
        print("Metric set to {0}".format(metric))

    if args.wt:
        if(args.wt == "no_weight"):
            global weight
            weight = None
            print("Weight set to {0}".format(weight))

k_nodes = None # default is None and hence all nodes are used
network_filename = "../HumanNetv2_converted.txt"
network_name = "HumanNetv2"
metric = "centrality"
weight = "weight"

def get_filename(network_name, k_nodes, type='pickle', metric='centrality', extension='pkl', weight=None):
    if k_nodes is None:
        filename_body = network_name
    else:
        filename_body = network_name+"_"+str(k_nodes)
    
    if weight is None:
        filename_body = filename_body +'_unweighted'
        
    if type == 'pickle':
        filename = "pickle_"+metric+"_dict_"+filename_body
    else:
        filename = filename_body+"_"+metric+"_values"
    
    filename = filename + "." + extension
    return filename
    
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

def calculate_metric_values(network_name, network_filename="", metric="centrality", k_nodes=None, weight='weight', seed=18):
    pickle_filename = get_filename(network_name, k_nodes, type="pickle", metric=metric, extension="pkl", weight=weight)
    print("Pickle file: {0}".format(pickle_filename))

    if os.path.exists(pickle_filename):
        f = open(pickle_filename, "rb")
        metric_values = pickle.load(f)
        print("Dictionary pickle loaded.")
    else:
        G = nx.read_weighted_edgelist(network_filename).to_undirected() # read undirected, weighted network
        if metric == "centrality":
            metric_values = nx.betweenness_centrality(G, k=k_nodes, normalized=True, weight=weight, 
                                                      endpoints=False, seed=seed)
        elif metric == "degree":
            metric_values = convert_list_of_tuples_to_dict(list(G.degree))
        
        print(metric_values)
        
        # pickle dictionary
        f = open(pickle_filename, "wb")
        pickle.dump(metric_values, f)
        print("Dictionary pickle created.")

    return metric_values

def convert_list_of_tuples_to_dict(list_of_tuples):
    dict = {}
    for item in list_of_tuples:
        dict[item[0]] = item[1]
    
    return dict
        
if __name__ == "__main__":
    parse_arguments()

    metric_values = calculate_metric_values(network_name=network_name, network_filename=network_filename, 
                                            weight=weight, metric=metric, k_nodes=k_nodes)

    results_csv_filename = get_filename(network_name, k_nodes, type='text', metric=metric, extension='csv', weight=weight)
    write_dict_to_CSV(metric_values, results_csv_filename, header=("gene,"+metric))
    
