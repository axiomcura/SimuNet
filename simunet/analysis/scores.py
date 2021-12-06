from collections import defaultdict
from itertools import permutations
from simunet.analysis.network import calculate_network_density
def _transform_to_list(data):
	"""Converts dictioanry to string

	Parameters
	----------
	data : dict
		data containing string and int as key value pair types

	Returns
	-------
	list
		list containing string repeated with their value

	Example
	-------
	>>> ds = {"net_1": 2, "net_2": 0, "net_3": 1}
	>>> ls = _dict_to_list(ds)
	>>> print(ls)
	>>> ["net_1", "net_1", "net_3"]
	"""
	# Type checking
	# -- keys must be string
	# -- values must be ints
	if not isinstance(data, dict):
		raise ValueError("data must be dict, {} was given".format(type(data)))
	elif not isinstance(list(data.keys()[0]), str):
		raise ValueError("keys must be a string, {} was given".format(type(list(data.keys()[0]))))
	elif not isinstance(list(data.values())[0], int):
		raise ValueError("values must be an int, {} was given".format(type(list(data.values()[0]))))

	ls = []
	for key, value in data.items():
		key_list = [key for i in range(value)]

		# update list
		ls += key_list
	return ls

def selection_score(mut_networks, string_db):
	"""Caculates selection score

	Parameters
	----------
	mut_networks : [type]
		[description]
	string_db :
	"""

	# creating a net and network score dictioanry
	net_scores = defaultdict(lambda: None)
	for idx, network in enumerate(mut_networks):
		key = "net_{}".format(idx+1)
		network_score = calculate_network_density(network)
		net_scores[key] = network_score

	ls_net_scores = _transform_to_list(net_scores)

	return ls_net_scores


def network_desnity_score(network, string_db):
	"""Caculates all edges in a network

	Parameters
	----------
	network : list
		gene network
	string_db : StringDB
		StringDB object

	Returns
	-------
	int
		number of edges in the whole network
	"""
	gene_pairs = permutations(network, 2)
	score = 0
	for gene1, gene2 in gene_pairs:
		paired = string_db.is_paired(gene1, gene2)
		if not paired:
			continue
		else:
			score += 1

	return score
