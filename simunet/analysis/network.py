import random
from simunet.analysis.methods import prix_fixe_selection, binned_prix_fixe_selction

def generate_fa_subnetworks(n, locus_data):
	"""Generating random subnetworks via prix-fixe method on
	FA gmt data

	Parameters
	----------
	n : int
		number of subnetworks generatred
	locus_data : dict
		gmt locus data

	Returns
	-------
	list
		array of arrays of different randomly generated subnetworks
	"""
	subnetworks = []
	for i in range(n):
		subnetwork = prix_fixe_selection(locus_data)
		subnetworks.append(subnetwork)
	return subnetworks


def generate_binned_subnetworks(n_networks, binned_data):
	"""generates random subnetworks using binned data

	Parameters
	----------
	n_networks : int
		number of subnetworks
	binned_data : dict
		bin id and genes as key value pairs

	Return
	------
	list
		populaitons of subnetworks
	"""
	random_binned_subnetworks = []
	for i in range(n_networks):
		network = binned_prix_fixe_selction(n_genes=12, binned_data=binned_data)
		random_binned_subnetworks.append(network)
	return random_binned_subnetworks


def generate_uninformative_pop(n_pop, n_networks, binned_genes):
	"""Generating n_populations with n_networks per populations. Uses
	binned genes for random selecting networks.

	Parameters
	----------
	n_pop : int
		number of populations
	n_networks : int, optional
		number of gene networks per populations
  		[by default 5000]

	Return
	------
	list
		list of n_populations that contains n_networks
	"""
	populations = []
	for pop_idx in range(n_pop):
		# generate randomly selected 5000 subnetworks using binned data -> returnsa list of list
		population = generate_binned_subnetworks(n_networks=n_networks, binned_data=binned_genes)
		populations.append(population)
		pass
	return populations


def mutation(locus_idx, gene, loci_data):
	"""Mutates genes to a different gene from the same locus

	Parameters
	----------
	locus_idx : [type]
		[description]
	gene : [type]
		[description]
	locus_data : dict
		locus index and associated genes key value pairs

	Returns
	-------
	str
		new gene
	"""
	locus_tag = "locus {}".format(locus_idx)
	genes_arr= loci_data[locus_tag]

	# preventing duplicate mutation
	while True:
		# if the randomly selected gene is the same as the input gene
		# -- pick again
		new_gene = random.choice(genes_arr)
		if new_gene == gene:
			new_gene = random.choice(genes_arr)
		else:
			break

	return new_gene

