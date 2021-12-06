# common.py
# used for basic file maniupations
import argparse
import random
from copy import deepcopy
from simunet.common.parsers import *
from simunet.analysis.methods import *


# TODO: Make sure that the returning list is not the same as the starting list
def simulate(subnetworks, locus_data , mut_freq=0.05):
	""" Applies genetic algorithm with a set of subnetworks at attempts
	to select select the best gene network.

	NOTE: read Tasan .et.al

	Returns
	-------
	list
		mutated subnetworks
	"""
	# NOTE: creating a copy for position updates
	# -- reassignment of python lists points to a refrence. Therefore any changes to
	# -- refrence will be applied to the original list
	mut_subnetworks = deepcopy(subnetworks)
	for network in mut_subnetworks:
		for locus_idx, gene in enumerate(network):
			if random.random() < mut_freq:
				print("mutation event on {} from locus {}".format(gene, locus_idx+1))
				new_gene = mutation(locus_idx+1, gene, locus_data)

				# index replacement
				# -- finding the index of the old gene and replace with new one
				network[network.index(gene)] = new_gene

	return mut_subnetworks


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


def generate_interaction_dict(subnetworks: list, stringDB: StringDB) -> dict:
	"""Generates a

	Parameters
	----------
	subnetworks : list
		list of generated subnet works
	stringDB : StringDB
		StringDB object

	Returns
	-------
	dict

	"""
	pass



if __name__ == '__main__':

    # CLI options
	# TODO: fix description
	# TODO: fix arguments designed for this program
	description = "simple program for generating SIF to full understand relationships between biological processes."
	parser = argparse.ArgumentParser(description=description)
	required = parser.add_argument_group("Required Arguments")
	required.add_argument("-i", "--input", type=str, metavar="INPUT",
						help="input file")
	parser.add_argument("-n", "--n_networks", type=int, required=False, default=5000,
						metavar="PARAMETER", help="generates n subnetworks")
	# parser.add_argument("-o", "--output", type=str, required=False, default="Simple_out",
	# 					metavar="PARAMETER", help="name of the outfile default='Simple_out'")
	# required.add_argument("-t", "--interaction_type", type=str, metavar="PARAMETER",
	# 					choices=["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"],
	# 					help="Interaction type")
	parser.add_argument('-db', "--database", type=str, metavar="FILE",
						help="path to database. Default path is `./Data/STRING.txt",
						default="./db/STRING.txt")
	args = parser.parse_args()

	# opening inputs
	db = StringDB(args.database)
	StringDB_df = db.to_pandas()
	loci_input = parse_input(args.input)

	# testing inputs
	# -- getting connecting genes and their counts
	# -- generating a proprocessed input (removing all 0 counts)
	connected_genes = get_connected_genes(StringDB_df)
	gene_counts = get_gene_counts(connected_genes)
	new_fa = preprocess_gmt(loci_input, gene_counts)

	# first step is go generated random networks
	print("MESSSAGE: Generating {} subnetworks ...".format(args.n_networks))
	subnetworks = generate_subnetworks(n=args.n_networks, locus_data=loci_input)

	print("Simulating with given subnetworks")
	mut_subnetworks = simulate(subnetworks, loci_input, mut_freq=0.5)

	print(subnetworks == mut_subnetworks)