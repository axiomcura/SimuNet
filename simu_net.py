# used for basic file maniupations
import argparse
import random
from copy import deepcopy
from statistics import mean

# import third parties
import pandas as pd

# simunet imports
from simunet.common.common import create_gmt_file
from simunet.common.parsers import *
from simunet.analysis.methods import *
from simunet.analysis.network import generate_fa_subnetworks, mutation, generate_uninformative_pop
from simunet.analysis.scores import *


def simulate(subnetworks, locus_data , string_db, mut_freq=0.05):
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
	if not isinstance(subnetworks, list):
		raise ValueError("Incorrect format, subnetworks should be in list format. You provided {}".format(type(subnetworks)))
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
	labled_networks = defaultdict(lambda: None)
	for network_idx, network in enumerate(mut_subnetworks):
		net_key = "net_{}".format(network_idx+1)
		for locus_idx, gene in enumerate(network):
			if random.random() < mut_freq:

				new_gene = mutation(locus_idx+1, gene, locus_data)

				# index replacement
				# -- finding the index of the old gene and replace with new one
				network[network.index(gene)] = new_gene
		labled_networks[net_key] = network
	# calculating selection scores
	selection_arr = selection_score_arr(mut_subnetworks, string_db)

	child_networks = []
	mate_counts = 0
	while True:

		if mate_counts == 5000:
			break

		# -- select random parent network from selection_arr
		# -- use the labled network to get their genes
		parent1 = random.choice(selection_arr)
		parent2 = random.choice(selection_arr)
		parent_1_network = labled_networks[parent1]
		parent_2_network = labled_networks[parent2]

		# -- using the mating function to get the child network
		# if parent!= parent_2:
		child_network = mate(parent_1_network, parent_2_network)
		child_networks.append(child_network)
		mate_counts += 1

	return child_networks


if __name__ == '__main__':

	# CLI options
	description = "Simunet: Program for optimizing random gene networks to find important gene interactions"
	parser = argparse.ArgumentParser(description=description)
	required = parser.add_argument_group("Required Arguments")
	required.add_argument("-i", "--input", type=str, metavar="INPUT",
						help="input file")
	parser.add_argument("-n", "--n_networks", type=int, required=False, default=5000,
						metavar="PARAMETER", help="generates n subnetworks, Default is 5000")
	parser.add_argument("-g", "--generations", type=int, required=False, default=30,
						metavar="PARAMETER", help="How many generations to simulate default is 30")
	parser.add_argument("-c", "--cutoff", type=float, required=False, default=0.05,
						 metavar="PARAMETER", help="difference in score cutoff. Ends program if covergence is detected. Default is 0.05")
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
	# -- bining data based on edge density
	connected_genes = get_connected_genes(StringDB_df)
	gene_counts = get_gene_counts(connected_genes)
	gene_ids = label_genes(gene_counts)
	binned_genes = bin_data(gene_counts, gene_ids)


	# first step is go generated random networks
	print("MESSSAGE: Generating {} subnetworks via prix fixe method...".format(args.n_networks))
	population = generate_fa_subnetworks(n=args.n_networks, locus_data=loci_input)
	org_population = deepcopy(population) # --> NOTE: used for score after?

	# generating uninformative networks
	print("MESSAGE: generating non-informative population")
	uninformative_pop = generate_uninformative_pop(n_pop=1000, binned_genes=binned_genes, n_networks=args.n_networks)

	print("MESSAGE: Optimizing FA associated randomly selected networks")
	# creating a for loop using  that controls the nmber of generaetions
	density_score = defaultdict(lambda: None)
	optmized_populations = []
	scores = []
	gen_pop_stats = []
	for gen_idx in range(args.generations):
		gen_key = "gen_{}".format(gen_idx+1)

		# Calculate statistics of the network
		stats = get_pop_stats(population, string_db=db)
		gen_pop_stats.append(stats)


		# creating next gen population and calulating its score
		next_gen_pop = simulate(population, loci_input, string_db=db, mut_freq=0.05)

		# calulate average network density
		avg_network_score = average_population_density_score(next_gen_pop, string_db=db)

		# checking for convergence
		# if lower that 0.05 change 3 conseq. times. break
		if len(scores) > 1:
			# print(gen_key, avg_network_score, abs(avg_network_score - scores[-1]))
			avg_change = abs(avg_network_score - scores[-1])
			if avg_change > 0.05:
				tracked_changes = []
			elif avg_change < 0.05:
				tracked_changes.append(avg_network_score)
				if len(tracked_changes) == 3:
					print("convergence at {} generations".format(gen_idx+1))
					break

		optmized_populations.append(next_gen_pop)
		scores.append(avg_network_score)
		population = next_gen_pop # -> updates the iterator with next gen pop

	print("MESSAGE: Calculating populations score of non-informative populations")
	pop_avg_score = defaultdict(lambda: None)
	for idx, population in enumerate(uninformative_pop):
		pop_key = "pop_{}".format(idx+1)
		network_scores = []
		for network in population:
			# calculate score
			network_score = calculate_gene_density_score(network=network, string_db=db)
			network_scores.append(network_score)
		pop_avg_score[pop_key] = np.mean(np.array(network_scores))


	print("MESSAGE: calculating p value")
	uniformative_scores = list(pop_avg_score.values())
	p_val = round(p_value_score(scores[-1], uniformative_scores), 4)
	print("p_value is {}".format(p_val))

	print("MESsAGE: Saving FA sta")
	columns = ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
	pd.DataFrame(gen_pop_stats, columns=columns).round(3).to_csv("GA_stats.txt", index=False)

	# creating gmt file
	print("MESSAGE: creating gmt file")
	create_gmt_file(loci_input, org_population, db)

	# top 10 networks of optomized population
	print("MESSAGE: saving top 10 networks of optomized population")
	opt_pop = optmized_populations[-1]
	opt_pop_scores = []
	for idx, network in enumerate(opt_pop):
		net_key = "Network{}".format(idx+1)
		score = calculate_gene_density_score(network, string_db=db)
		result = [net_key, score]
		opt_pop_scores.append(result)

	cols = ["network", "score"]
	top_10 = pd.DataFrame(opt_pop_scores, columns=cols).sort_values(by=["score"], ascending=False).iloc[:10]
	print(top_10)
	top_network = top_10.iloc[0]["network"]
	outname = 'Day3_Output_{}_pval{}.txt'.format(top_network, p_val)
	top_10.to_csv(outname, index=False)

	print("complete!")
