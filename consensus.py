#!/usr/local/bin/python
"""
Module to generate consensus sequence
"""
from poaligner import *
import sys

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

	return node['previous']


def score_assignment(dag):
	"""
	Initializes node scores and edge weights
	Computes & assigns best score for all nodes in DAG 
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


def get_max_score_node(dag):
	"""
	return the node with maximum score
	"""
	max_node = dag[0]
	for node in dag[1:]:
		max_node = node if node['score'] > max_node['score'] else max_node

	return max_node


def do_consensus(dag):
	"""
	returns consensus sequence using back-pointers from node with max score
	"""
	cons_seq = ""
	max_node = get_max_score_node(dag) # have to try alternatives
	cons_seq += max_node['character']
	while max_node['previous']:
		prev_node = max_node['previous']
		cons_seq += dag[prev_node]['character']
		max_node = dag[prev_node]
	cons_seq += max_node['character']

	return cons_seq[::-1]

if __name__ == '__main__':
	if len(sys.argv) < 2:
		exit("Required format: python consensus.py <.po file>")
	po_msa_file = sys.argv[1]
	dag = convert_po_msa_to_dag(po_msa_file)

	# assigning edge-weights & node scores
	score_assignment(dag)
	
	# generating consensus
	sequence = do_consensus(dag)
	print (sequence)
