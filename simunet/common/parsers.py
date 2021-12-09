import sys
from collections import defaultdict
import pandas as pd

class StringDB:
    """ Creates a simple python object that handles lookups into the
    STRING.txt file.
    Arguments
    ---------
    fname : str
        path to database
    Returns
    -------
    StringDB Object
        Small container where the the String.txt is stored
    """

    def __init__(self, fname):
        db = self._parse_string_db(fname)
        self._db = db
        self.fname = fname
        self.size = "{} MB".format(round(sys.getsizeof(db)/1024**2, 4))


    def _parse_string_db(self, fname: str) -> dict:
        """ Parses the STRING.txt file and converts it into a dictionary
        Arguments
        fname : str
            path to STRING.txt file
        Returns
        -------
        dict
            dictionary containing protein pairs as key and scores as its value
        # Example:
        >>> results = {"PROT1 PROT2" : score, ...}
        """

        # Reads the STRING.txt file and iterates over all lines
        # each line is split into two parts
        # -= protein pairs (key)
        # -- score (value)
        # then stored into dictionary
        db_contents = defaultdict(lambda: None)
        with open(fname, "r") as infile:
            record_contents = infile.readlines()
            for records in record_contents:
                data = records.replace("\n", "").split("\t")
                protein_pair = "{} {}".format(data[0], data[1])
                score = data[2]
                db_contents[protein_pair] = score
        return db_contents

    def to_pandas(self) -> pd.DataFrame:
        """ Converts the StringDB object into a pandas object"""
        return pd.read_csv(self.fname, delimiter="\t", names=["gene1", "gene2", "score"])

    def get_paired_score(self, gene1, gene2):
        """Interaction score between two genes. If query is not found,
        a score of 0 is returned"

        Parameters
        ----------
        gene1 : str
            first gene for query
        gene2 : str
            second gene for query

        Return
        ------
        float
            gene density score in network
        """
        query = "{} {}".format(gene1, gene2)
        score = self._db[query]
        if score is None:
            return 0.0
        return float(score)


    def is_paired(self, gene1, gene2):
        """ checks if the two genes are paired

        Paramters
        ---------
        gene1 : str
            query of gene 1
        gene2 : str
            query of gene 2

        Returns
        -------
        bool
            Returns True if there's a connection, False if no connection exists
        """
        query = "{} {}".format(gene1, gene2)
        check = self._db[query]
        if check is None:
            return False
        return True


def preprocess_gmt(fa_input: dict, gene_counts: dict) -> dict:
	"""Removes genes from each locus where its count is zero when searched in the StringDB

	Parameters
	----------
	fa_input : dict
		original fa_input
	gene_counts : dict
		dictionary containing gene counts

	Results
	-------
	dict
		preprocesed locus and gene array as key value pairs
	"""
	prep_fa = defaultdict(lambda: None)
	for idx, (locus_name, gene_arr) in enumerate(fa_input.items()):
		locus_name = "locus {}".format(idx+1)
		gene_list = []
		for gene in gene_arr:
			try:
				count = gene_counts[gene]
			except KeyError:
				continue
			if count == 0:
				continue
			gene_list.append(gene)
		prep_fa[locus_name] = gene_list

	return prep_fa


def parse_input(input_file : str) -> dict:
    """ Parses input file and converts it into a dictionary
    Arguments
    ---------
    input_file : str
        path to input file
   Returns
   -------
   dict
        dictionary containing the locus name as the key and all
        the genes within the locus as value

   Example
   --------
    >>> # this is an example output of this function.
    >>> parsed_input = {"locus 1" : ["gene1", "gene2", "gene3", "geneN"],
    >>>                 "locus 2" : ["gene1", "gene2", "gene3", "geneN"]}
    """

    locus_genes = defaultdict(lambda: None)
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for locus_idx, line in enumerate(lines):
            data = line.strip("\n").split("\t")
            locus = "locus {}".format(locus_idx+1)
            genes = data[2:]
            locus_genes[locus] = genes

    return locus_genes


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