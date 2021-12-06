import random
from collections import defaultdict
import pandas as pd

def prix_fixe_selection(locus_data):
	"""Selects randomly genes within each locus and returns a list of genes
	that represents the one sub network

	Parameters
	----------
	locus_data : dict
		contains locus and associated genes data. Where a row represents one locus
		containing the associated genes

	Returns:
	--------
	list
		containing a list of genes
	"""
	subnetwork = []
	for gene_list in locus_data.values():
		random_gene = random.choice(gene_list)
		subnetwork.append(random_gene)
	return subnetwork




def generate_binned_subnetworks(n, nbins):
    """Generating random subnetworks via prix-fixe method on
    """
    pass

def get_connected_genes(string_df: pd.DataFrame):
	"""
	Searches within STRING database and obtains informations genetic connectivity.

	Arguments:
	---------
	string_df : pd.DataFrame
		Pandas DataFrame of the STRING database

	Returns
	-------
	dict
		Gene name and array of associated genes as key value pairs
	"""
	# generating a set of all known genes within string database
	# -- both columns genes together and converts it into a set (removes all repeats)

	# grouping data based on gene name
	grouped_df = string_df.groupby(by=["gene1"])
	connected_genes = defaultdict(None)
	for gene, gene_df in grouped_df:
		conn_genes_arr = gene_df["gene2"].values.tolist()
		connected_genes[gene] = conn_genes_arr

	return connected_genes


def get_gene_counts(connected_genes_data):
	"""
	Creates a sorted dictioanry of genes and its edge density.

	Arguments
	---------
	connected_genes_data : dict
		Dictioanry containing genes name as the key and all the other associated genes as the values

	Returns
	-------
	dict
		Gene name and edge desnity key value pairs,
	"""
	gene_counts = defaultdict(None)
	for gene, connected_genes in connected_genes_data.items():
		gene_counts[gene] = len(connected_genes)

	# sorts the genes by edge density
	sorted_gene_counts = {k: v for k, v in sorted(gene_counts.items(), key=lambda item: item[1])}

	return sorted_gene_counts


def label_genes(counts):
    """ Provides genes with a unique number as and id

    Arguments
    ---------
	counts : dict
		Gene counts dictionary where each gene name has as asscoiated value, which is edge density

	Returns
	--------
	dict
		labled genes with unique ids
    """
    labeled_genes = defaultdict(lambda: None)
    for idx, gene_name in enumerate(counts.keys()):
        labeled_genes[idx] = gene_name
    return labeled_genes

