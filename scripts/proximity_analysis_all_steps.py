import networkx as nx
import pandas as pd
import time
import random
import numpy
import datetime
from multiprocessing import Process

print("Script started...")

# ============================== FUNCTIONS ==============================
def has_path(G,node_from,node_to): # returns true if the graph G has a path from source node(node_from) to target node (node_to)
    return(nx.has_path(G,node_from,node_to))


def get_shortest_path_length_between(G, source_id, target_id): #calculates shortest path length between source and target nodes in the graph G
    return nx.shortest_path_length(G, source_id, target_id) 

'''
The function below uses the previous two functions to calculate the path length 
between each source node (drug target) and its closest disease target.
See Guney 2015 for more details on the approach (Network-based in silico drug efficacy screening).
'''
def calculate_closest_distance(network, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from: #nodes_from is a list of drug targets (eg. from drugbank, chembl, etc.)
        values = [] # will store the shortest path length between a source node and all disease targets here
        for node_to in nodes_to: #nodes_to is a list of known disease targets
            # print("from - to", node_from, node_to)
            if not has_path(network,node_from,node_to): continue
            val = get_shortest_path_length_between(network, node_from, node_to)
            values.append(val)
        if len(values) == 0:    continue
        d = min(values) # the shortest path between a source node and its closest disease target
        # print (d)
        values_outer.append(d)
    closest_d = numpy.mean(values_outer) # the average shortest path length between any source node (drug target) and its closest disease target
    # print (d)
    return closest_d


def get_degree_binning(g, bin_size):
    '''
    This function creates the bins from a network g.
    Starting from a list of nodes with the lowest degree, it adds nodes with the same degree to the bin until it reaches the set bin size.
    If number of nodes with some degree is lower then bin size, it combines with other nodes with degree + 1 to meet bin size.
    '''
    degree_to_nodes = {}
    # the below two lines compute the degree of each node in the graph.
    # the setdefault method is used to add each node to a list of nodes with the same degree in the dictionary
    for node, degree in g.degree(): 
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values = sorted(values) # values becomes a list, sorted from lowest to highest degree
    bins = []
    i = 0 # this is the iterator that iterates over each degree, starting from the first item in the list (lowest degree)
    while i < len(values):
        low = values[i] # this is the i-th degree in the values list
        val = degree_to_nodes[values[i]] # a list of the nodes with i-th degree (low)
        while len(val) < bin_size:
            # while the number of nodes in a bin is lower than the bin size, than nodes with degree i+1 will be added to the bin
            # bin size is chosen by the user - in the paper this is set to 100
            i += 1 # next iteration (move to the next degree in the list)
            if i == len(values): # breaks when the last item in the list is reached
                break
            # starting from a list of nodes with the lowest degree, it adds nodes with degree lowest + 1 to the val list until it reaches the set bin size
            val.extend(degree_to_nodes[values[i]]) # val will be extended with the next set of nodes with degree i+1 (low +1).
        if i == len(values):
            i -= 1
        high = values[i] # this is the highest degree
        i += 1
        # print i, low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins

# this function lists all nodes in the same bin as each seed node
def get_degree_equivalents(seeds, bins, g): 
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed) #extract degree of the seed node
        for l, h, nodes in bins: #it takes low, high degree and nodes for each bin
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

    
def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False,
                                        seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes() # list of nodes in the network
    for i in range(n_random): # decided by the user (how many times will the random iterations be repeated?) usually this is = 1000
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            # the lines below pick random nodes matching the degree (same bin) of the real nodes
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network) # lists nodes in the same bin as the node of interest
            # now choose a random node from the same bin as the real node
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times) - to make sure it doesn't choose the same node
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    '''
    This function creates a n_random number of lists of random nodes with the same degree binning as the real nodes (when degree_aware=True).
    usually n_random = 1000 because we often do 1000 iterations.
    '''
    if bins is None:
        # Get degree bins of the network (if they aren't already supplied
        bins = get_degree_binning(network, min_bin_size)
    # pick the random nodes
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware,
                                                                         seed=seed)
    return nodes_random

def calculate_proximity(network, drug, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None,
                        n_random=1000, min_bin_size=100, seed=452456):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    """

    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network # select only nodes_from (drug targets) that are located in the network
    nodes_to = set(nodes_to) & nodes_network # select only nodes_to (disease targets) that are located in the network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None  # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to) # this is the real distance
    
    # now do 1000 iterations using random nodes
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size)
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins=bins, n_random=n_random,
                                             min_bin_size=min_bin_size, seed=seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins=bins, n_random=n_random, min_bin_size=min_bin_size,
                                           seed=seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random))  # n_random
    # now calculates the closest distance using random nodes. Repeat x1000
    for i, values_random in enumerate(random_values_list):
        #print('iteration ', i)
        nodes_from, nodes_to = values_random
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to)
    m, s = numpy.mean(values), numpy.std(values) # do mean and stdev of random iterations
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    dict = {'drug': drug, 'distance': d, 'z_score': z}
    return dict


# ============================== PROXIMITY ANALYSIS ==============================
# Load PPI
print("Loading PPI network...")
ppi_df = pd.read_csv("../data/networks/combined_PPI.csv", sep=",")
ppi_graph = nx.from_pandas_edgelist(ppi_df, "GeneA", "GeneB")

# Load DPI
print("Loading DPI...")
dpi_df = pd.read_csv("../data/networks/combined_DPI.csv", sep=",")
dpi_df.dropna(inplace=True)
dpi_df = dpi_df[dpi_df["Drug_Target"].isin(ppi_graph.nodes)] # Keep only drug targets that are in the PPI graph
drug_list = dpi_df["Drug_Name"].unique()
print("Preparing drug list and bins...")
drug_to_targets = dpi_df.groupby("Drug_Name")["Drug_Target"].apply(set).to_dict() # This is a mapping of drug names to their targets

# Get degree bins
bins = get_degree_binning(ppi_graph, bin_size=100)

print("Starting multiprocessing...")

# Proximity function
def run_proximity_analysis(step_name, deg_file_path, output_csv_path, ppi_graph, drug_to_targets,
                           bins, n_random=1000, min_bin_size=100, seed=452456):
    """
    Calculate proximity for step_name using drugs and targets in drug_to_targets and
    disease genes in deg_file_path
    Saves the result as a CSV file at output_csv_path
    """
    # Open a step-specific log file
    log_path = f"{step_name}_log.txt"
    def log(message):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(log_path, "a") as f:
            f.write(f"[{timestamp}] {message}\n")

    log(f"===== Processing {step_name} =====")
    
    # ========== LOAD DISEASE GENES ==========
    disease_genes = pd.read_csv(deg_file_path, header=None)
    log(f"{len(disease_genes)} disease genes in {step_name}")
    disease_genes = disease_genes.iloc[:,0].dropna().unique().tolist()
    # Keep only disease genes that are in the PPI graph
    disease_genes = [gene for gene in disease_genes if gene in ppi_graph.nodes]
    log(f"{len(disease_genes)} disease genes in {step_name} after filtering for PPI graph")

    # ========== GENERATE RANDOM SETS OF DISEASE GENES ==========
    # For each disease gene, pick a random node from the same degree bin, and repeat this 1000 times,
    # generating 1000 sets of random disease genes
    random_disease_genes = get_random_nodes(
        disease_genes,
        ppi_graph,
        bins=bins,
        n_random=n_random,
        min_bin_size=min_bin_size,
        seed=seed,
        degree_aware=True)
    log(f"Generated {len(random_disease_genes)} sets of random disease genes")

    # ========== RUN PROXIMITY ANALYSIS ==========
    results = []
    # Iterate over each drug and calculate proximity to disease genes
    for i, (drug, targets) in enumerate(drug_to_targets.items()):
        start = time.time()
        log(f"Processing drug: {drug}")

        # Calculate proximity to disease genes
        result = calculate_proximity(
            ppi_graph,
            drug,
            targets,
            disease_genes,
            nodes_from_random=random_disease_genes,
            bins=bins,
            n_random=n_random,
            min_bin_size=min_bin_size,
            seed=seed
        )

        if result is not None:
            results.append(result)
        else:
            log(f"No proximity data for drug: {drug}")

        if i % 50 == 0:
            log(f"{i}/{len(drug_to_targets)} drugs processed")
        
        end = time.time()
        log(f"Runtime for {drug}: {end - start:.2f} sec")

    # Save and display
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv_path, index=False)
    log(f"Saved results to {output_csv_path}")
    log(f"Top 10 prioritised drugs:\n{results_df.sort_values('z_score').head(10)}")

# Define DEG file paths and output paths
steps = {
    "step1": {
        "deg_file": "../data/networks/step1_deg.csv",
        "output_csv": "../results/humanPVATsn/network_analysis/proximity_step1.csv"
    },
    "step2": {
        "deg_file": "../data/networks/step2_deg.csv",
        "output_csv": "../results/humanPVATsn/network_analysis/proximity_step2.csv"
    },
    "step3": {
        "deg_file": "../data/networks/step3_deg.csv",
        "output_csv": "../results/humanPVATsn/network_analysis/proximity_step3.csv"
    },
    "full": {
        "deg_file": "../data/networks/full_diff_deg.csv",
        "output_csv": "../results/humanPVATsn/network_analysis/proximity_full.csv"
    }
}

if __name__ == "__main__":
    processes = []
    for step, path in steps.items():
        p = Process(target=run_proximity_analysis, kwargs=dict(
            step_name=step,
            deg_file_path=path["deg_file"],
            output_csv_path=path["output_csv"],
            ppi_graph=ppi_graph,
            drug_to_targets=drug_to_targets,
            bins=bins,
            n_random=1000,
            min_bin_size=100,
            seed=452456
        ))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()