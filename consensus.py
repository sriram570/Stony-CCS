#!/usr/local/bin/python
"""
Module to generate consensus sequence
"""
from poaligner import *
import random
import sys

# Allowed scoring functions and traversal algorithms
SCORING_FUNCTIONS = ["edge_weight_based_score",
					 "pb_like_score"]
TRAVERSAL_ALGOS   = ["max_score",
					 "max_in_edges",
					 "max_sequences",
					 "random",
					 "optimal_random"]

# ========================================== INITIALIZATIONS =================================================

def init_score(dag):
	"""
	Initializing scores of nodes in DAG
	"""
	for node in dag:
		node['score'] = 0


def init_edge_weights(dag):
	"""
	Assigning weights to edges in DAG
	Returns a dict with keys as (start_node_index, end_node_index) and 
								values as edge-weights b/w these indexes
	"""
	edge_weights = {}
	for node in dag:
		for prev_node in node['incoming']:
			edge_weights[str(prev_node), str(node['index'])] = \
				len(dag[prev_node]['sequences'].intersection(node['sequences']))

	return edge_weights


# ========================================= SCORING FUNCTION 1 ===============================================

def best_predecessor_node(dag, node, incoming_nodes, edge_weights):
	"""
	returns index of the best predecessor node
	"""
	max_weight = 0
	for neighbor in incoming_nodes:
		if edge_weights[str(neighbor), str(node['index'])] > max_weight:
			max_weight = edge_weights[str(neighbor), str(node['index'])]
			node['previous'] = neighbor

		# when edge weights are same, takes predecessor scores into consideration (can try other alternatives)
		elif edge_weights[str(neighbor), str(node['index'])] == max_weight:
			if (dag[node['previous']]['score'] + edge_weights[str(node['previous']), str(node['index'])]) < \
					(dag[neighbor]['score'] + edge_weights[str(neighbor), str(node['index'])]):
				node['previous'] = neighbor

	if len(incoming_nodes) == 0:
		node['previous'] = None

	return node['previous']


def score_assignment(dag):
	"""
	scoring function to initialize node scores and edge weights
	Computes & assigns best score for all nodes in DAG
	scoring ftn --> score(node) = score(best_prev_node) + edge_weight(best_prev_node, curr_node)  
	"""
	init_score(dag)
	edge_weights = init_edge_weights(dag)

	# assigning start node's score
	dag[0]['score'] = 0

	# computing scores of other nodes
	for node in dag[1:]:
		index = node['index']
		if len(node['incoming']) > 0:
			best_prev_node = best_predecessor_node(dag, node, node['incoming'], edge_weights)
			node['score'] = dag[best_prev_node]['score'] + edge_weights[str(best_prev_node), str(node['index'])]
		elif len(node['incoming']) == 0:
			node['previous'] = None

# ========================================= SCORING FUNCTION 2 ================================================

def set_out_edges(dag):
	"""
	this method computes and sets the #outgoing edges for each node in POA graph
	"""
	for node in dag:
		node['outgoing'] = 0

	for node in dag:
		for prev_node in node['incoming']:
			dag[prev_node]['outgoing'] += 1


def pb_like_score_assignment(dag):
	"""
	scoring function similar to PacBio's scoring function
	scoring ftn --> score(node) = 2*contained_reads - #outgoing_edges
	"""
	init_score(dag)
	edge_weights = init_edge_weights(dag)
	set_out_edges(dag)

	for node in dag:
		node['score'] = 2*len(node['sequences']) - node['outgoing']
		best_predecessor_node(dag, node, node['incoming'], edge_weights) # computes back-pointer


# ========================================= INVOKING SCORING FUNCTIONS =========================================

def scoring_function(dag, scoring_func):
	if scoring_func == "edge_weight_based_score":
		score_assignment(dag)
	elif scoring_func == "pb_like_score":
		pb_like_score_assignment(dag)


# ================================ MULTIPLE TRAVERALS TO OBTAIN CONSENSUS  ====================================

def get_max_score_node(dag):
	"""
	returns the node with maximum score
	"""
	max_node = dag[0]
	for node in dag[1:]:
		max_node = node if node['score'] >= max_node['score'] else max_node

	return max_node


def max_in_edges_node(dag):
	"""
	returns the last node with maximum incoming edges in the topologically sorted POA graph
	"""
	best_node = dag[0]
	for node in dag[1:]:
		best_node = node if len(node['incoming']) >= len(best_node['incoming']) else best_node

	return best_node


def max_seq_node(dag):
	"""
	returns the last node containing maximum reads in the topologically sorted POA graph
	"""
	best_node = dag[0]
	for node in dag[1:]:
		best_node = node if len(node['sequences']) >= len(best_node['sequences']) else best_node

	return best_node


def random_node(dag):
	"""
	returns any node having score in range [max_score - threshold, max_score]
	"""
	node_list = []
	max_node = get_max_score_node(dag)
	lower_range, upper_range = max_node['score'] - 0.01*max_node['score'], max_node['score']

	# finds a subset of potential end nodes
	for node in dag:
		if node['score'] >= lower_range and node['score'] <= upper_range:
			node_list.append(node)

	return random.choice(node_list)


def optimal_random_node(dag):
	"""
	extracts a subset of nodes having score in range [max_score - threshold, max_score]
	from the extracted subset, returns the last node containing maximum reads
	"""
	node_list = [] # subset of nodes
	max_node = get_max_score_node(dag)
	lower_range, upper_range = max_node['score'] - 0.01*max_node['score'], max_node['score']

	# finds a subset of potential end nodes
	for node in dag:
		if node['score'] >= lower_range and node['score'] <= upper_range:
			node_list.append(node)

	# from the filtered subset, selects an optimal node containing max reads
	opt_node = node_list[0]
	for node in node_list[1:]:
		if len(node['sequences']) >= len(opt_node['sequences']):
			opt_node = node

	return opt_node


# ========================================= CONSENSUS GENERATION ==============================================

def do_consensus(dag, traversal_algo):
	"""
	returns consensus sequence using back-pointers from a selected node
	"""
	cons_seq = ""

	# node selection alternatives for POA traversal
	if traversal_algo == "max_score":
		max_node = get_max_score_node(dag)
	elif traversal_algo == "max_in_edges":
		max_node = max_in_edges_node(dag)
	elif traversal_algo == "max_sequences":
		max_node = max_seq_node(dag)
	elif traversal_algo == "random":
		max_node = random_node(dag)
	elif traversal_algo == "optimal_random":
		max_node = optimal_random_node(dag)

	cons_seq += max_node['character']
	while max_node['previous']:
		prev_node = max_node['previous']
		#print max_node['previous']
		cons_seq += dag[prev_node]['character']
		max_node = dag[prev_node]
	cons_seq += max_node['character']

	return cons_seq[::-1]


def consensus_to_fasta(sequence, out_file=None):
	"""
	write consensus to fasta file "consensus.fa" for blasr comparison
	"""
	f = open(out_file, 'w')
	f.write(">consensus\n")
	f.write(sequence)
	f.close()


if __name__ == '__main__':
	if len(sys.argv) < 2:
		exit("Required format: python consensus.py <.po file>")
	po_msa_file = sys.argv[1]
	dag = convert_po_msa_to_dag(po_msa_file)

	# assigning edge-weights & node scores
	scoring_function(dag, "pb_like_score")

	# generating consensus
	sequence = do_consensus(dag, "max_in_edges")
	consensus_to_fasta(sequence, "consensus.fa")