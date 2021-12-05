# common.py
# used for basic file maniupations
import argparse
from simunet.common.parsers import *
from simunet.analysis.methods import *


def simulate():
	""" Applies genetic algorithm with a set of subnetworks at attempts
	to select select the best gene network.

	NOTE: read Tasan .et.al
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
	fa_input = parse_input(args.input)

	# testing inputs
	# -- getting connecting genes and their counts
	# -- generating a proprocessed input (removing all 0 counts)
	connected_genes = get_connected_genes(StringDB_df)
	gene_counts = get_gene_counts(connected_genes)
	new_fa = preprocess_gmt(fa_input, gene_counts)

	# first step is go generated random networks
	print("MESSSAGE: Generating {} subnetworks ...".format(args.n_networks))
	subnetworks = generate_subnetworks(n=args.n_networks, locus_data=fa_input)

