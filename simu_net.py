# common.py
# used for basic file maniupations
import argparse
import random
from copy import deepcopy
from simunet.common.parsers import *
from simunet.analysis.methods import *
from simunet.analysis.network import generate_subnetworks, mutation
from simunet.analysis.scores import selection_score_arr, average_population_density_score

def simulate(subnetworks, locus_data , string_db, generations=None, cutoff=None, mut_freq=0.05):
	"""Applies genetic algorithm with a set of subnetworks at attempts
	to select  the best gene network. Convergence is captured after
	+3 generation.

	Parameters
	----------
	subnetworks1 : list
		subnetworks within the first population
	locus_data : dict
		locus name and associated genes as key value pairs
	string_db : dict
		gene pairs and their respective scores as key value
		pairs
	generations : int
		how many generations should the simulation run
	cutoff : float
		cut off value when detecting convergence
	mut_freq : float, optional
		mutation rate (default 0.05)
	"""
	# type checking
	if not isinstance(subnetworks, list) and not isinstance(subnetworks, list):
		raise ValueError("Incorrect format, subnetworks1 and subnetworks2 shold be lists. You prvided {}, {}".format(type(subnetworks1), type(subnetworks2)))
	elif not isinstance(locus_data, dict) and not isinstance(string_db, dict):
		raise ValueError("Incorrect format, locus_data and string_db shold be dicts. You prvided {}, {}".format(type(locus_data), type(string_db)))
	elif mut_freq > 1.0:
		raise ValueError("mutation frequences my be lower than 1.0, you provided {}".format(mut_freq))

	# NOTE: creating a copy for position updates
	# -- reassignment of python lists points to a refrence. Therefore any changes to
	# -- refrence will be applied to the original list
	avg_density_score = []
	result_child_networks = defaultdict(lambda: None)

	mut_subnetworks = deepcopy(subnetworks)
	for network1 in mut_subnetworks[:10]:
		for locus_idx, gene in enumerate(network1):
			if random.random() < mut_freq:
				print("mutation event on {} from locus {}".format(gene, locus_idx+1))
				new_gene = mutation(locus_idx+1, gene, locus_data)

				# index replacement
				# -- finding the index of the old gene and replace with new one
				network1[network1.index(gene)] = new_gene


	# calculating selection scores
	# TODO: test algorithm
	print("Testing selection ")
	selection_arr = selection_score_arr(mut_subnetworks, string_db)

	child_networks = []
	mate_counts = 1
	while True:
		if mate_counts == 5000:
			break
		# -- select random parent network from selection_arr
		parent_1 = random.choice(selection_arr)
		parent_2 = random.choice(selection_arr)

		# -- using the mating function to get the child network
		if parent_1 != parent_2:
			child_network = mate(parent_1, parent_2)
			child_networks.append(child_network)
			mate_counts += 1


	return child_networks


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
	# parser.add_argument("-g", "--generations", type=int, required=False, default=20,
	# 					metavar="PARAMETER", help="How many generations to simulate")
	# parser.add_argument("-c", "--cutoff", type=float, required=False, default=0.1,
    #                  	metavar="PARAMETER", help="difference in score cutoff. Ends program if covergence is detected")
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
	population = generate_subnetworks(n=args.n_networks, locus_data=loci_input)
	org_population = deepcopy(population) # --> NOTE: used for score after?

	print("Simulating with given subnetworks")
	# creating a for loop using  that controls the nmber of generaetions
	density_score = defaultdict(lambda: None)
	populations = defaultdict(lambda: None)
	population = None

	# TODO: change that value within range function into args.generations argument
	for gen_idx in range(10):
		gen_key = "gen_{}".format(gen_idx)

		# creating next gen population and calulating its score
		next_gen_pop = simulate(population, loci_input, string_db=db, mut_freq=0.5)

		# calulate average network density
		avg_network_score = average_population_density_score(next_gen_pop, string_db=db)
		density_score[gen_key] = avg_network_score
		populations[gen_key] = next_gen_pop
		population = next_gen_pop # -> updates the iterator with next gen pop