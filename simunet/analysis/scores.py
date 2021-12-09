from collections import defaultdict
from itertools import combinations


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
	elif not isinstance(list(data.keys())[0], str):
		raise ValueError("keys must be a string, {} was given".format(type(list(data.keys()[0]))))
	elif not isinstance(list(data.values())[0], int):
		raise ValueError("values must be an int, {} was given".format(type(list(data.values()[0]))))

	ls = []
	for key, value in data.items():
		if value == 0:
			continue
		key_list = [key for i in range(value)]

		# update list
		ls += key_list
	return ls


def selection_score_arr(network, string_db):
	"""Caculates selection score

	Parameters
	----------
	network : list
		list of genes with in a network
	string_db : dict
		gene pairs and scores as key value pairs

	Returns
	-------
	list
		array of network names that is repeated.
	"""

	# creating a net and network score dictioanry
	net_scores = defaultdict(lambda: None)
	for idx, network in enumerate(network):
		key = "net_{}".format(idx+1)
		network_score = network_edge_desnity_score(network, string_db)
		net_scores[key] = network_score**3

	ls_net_scores = _transform_to_list(net_scores)

	return ls_net_scores


def network_edge_desnity_score(network, string_db):
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
	gene_pairs = combinations(network, 2)
	matches = []
	for gene1, gene2 in gene_pairs:
		is_paired = string_db.is_paired(gene1, gene2)
		if not is_paired:
			continue
		matches.append(is_paired)

	score = len(matches)
	return score

def calculate_gene_density_score(network, string_db):
	"""Calculates the gene score within each network

	Parameters
	----------
	network : list
		list of genes in a network
	string_db : stringDB
		stringDB object

	Returns
	-------
	int
		Sum of gene scores within a network
	"""

	scores = []
	gene_pairs = combinations(network, 2)
	for gene1, gene2, in gene_pairs:
		score = string_db.get_paired_score(gene1, gene2)
		scores.append(score)
	total_score = sum(scores)
	return total_score


def average_population_density_score(population, string_db):
	"""calculates the average population density score by
 	calcuating the gene weights per each network


	Parameters
	----------
	population : list
		list of list containing networks
	string_db : [type]
		[description]
	"""

	population_density_scores = []
	for network in population:
		network_score = calculate_gene_density_score(network, string_db)
		population_density_scores.append(network_score)

	avg_population_score = sum(population_density_scores)/len(population_density_scores)
	return avg_population_score



