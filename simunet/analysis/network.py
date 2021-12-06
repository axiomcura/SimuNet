import random
from simunet.analysis.methods import prix_fixe_selection

def generate_subnetworks(n, locus_data):
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

def mate():
	"""[summary]
	"""
	pass
