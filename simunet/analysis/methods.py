import random
from collections import Counter
from collections import defaultdict
import pandas as pd

from simunet.common.errors import NetWorkError

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


def mate(parent1, parent2):

    """[summary]

    Parameters
    ----------
    parent1 : [type]
        list of genes associated
    parent2 : list
        list of genes associated with parent1
    """
    # randomly selecting 6 genes from
    parent1_genes = random.choices(parent1, k=6)
    parent2_genes = random.choices(parent2, k=6)
    child_network = parent1_genes + parent2_genes

    return child_network
    # TODO: should I remove duplicated genes from child networks? 
    # checking for duplicate genes
    # NOTE: may not be necessary
    # n_dups = _check_duplicates(child_network)
    # if n_dups is not None:
    #     print("WARNING: {} Duplicated genes found in resulting child subnetwork".format(n_dups))
    #     child_network = _remove_duplicates(child_network)


def _check_duplicates(child_network):
    """ Removes duplicated genes in a subnetwork"""
    if len(child_network) == len(set(child_network)):
        return None
    else:
        n_dups =  (len(child_network) - len(set(child_network))) / 2
        return n_dups


def _remove_duplicates(child_network, parent1_network, parent2_network):
    print("MESSAGE: resampling genes...")
    # continously removing duplicates until a unique network is found
    attempts = 0
    while True:
        if attempts == 100:
            raise NetWorkError("Unable to remove duplicate genes from resulting child network")

        # find the indexes of these duplicated genes
        gene_name, counts = Counter(child_network).most_common()[0]
        dup_index = [i for i,value in enumerate(child_network) if value == gene_name]
        for idx in dup_index:
            if dup_index < 6:
                # randomly select from parent 2
                new_gene = random.choices(parent1_network)
                child_network[idx] = new_gene
            elif dup_index > 6:
                # randomly select gene from parent 2
                new_gene = random.choices(parent2_network)
                child_network[idx] = new_gene

    # research for duplicates genes in child network use counter
        # if max value is > 2,
            # increment attempt + 1
            # resample again

        # return child list if not duplicates found


