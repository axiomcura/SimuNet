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


    def get_pair_score(self, locus_name: str, proteins: list, interaction_type="pp") -> list:
        """ Accepts a list of genes and queries to the database
        Summary:
        -------
        Creates all possible permutations of protein interactions within the protein list.
        Each pair will be queried into database and returns a score. The database will return
        'None' if the score is not found and will not be recorded into the results. In addition,
        a "tracking" list is also implemented to prevent repetitive query. This means that
        reversed interaction queries will be ignored if the original query has been recorded.
        Argument
        -------
        genes : list
            list of genes found in the locus
        interaction_type : str (choices=["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"])
            Interaction type of both genes/molecules. Supported interaction types are:
            - p:protein - protein interaction
            - pd: protein -> DNA
            - pr: protein -> reaction
            - rc: reaction -> compound
            - cr: compound -> reaction
            - gl: genetic lethal relationship
            - pm: protein-metabolite interaction
            - mp: metabolite-protein interaction
        More information about compatibility found here:
        http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html#sif-format
        Returns
        -------
        list
            Contains a list of strings that describes the interaction type
            between two genes and its score. This cotnents is what is going
            to be used to produce the sif file
        """

        # type checking
        if not isinstance(proteins, list):
            genes = [proteins]

        # getting all possible combinations
        # -- query all combinations into the database
        # -- reversed queries are ignored if original query is recorded
        pairs = permutations(proteins, 2)

        results = []
        searched = []
        for gene1, gene2 in pairs:
            query = "{} {}".format(gene1, gene2)
            query_rev = "{} {}".format(gene2, gene1)
            score = self._db[query]
            if score == None:
                continue
            if  query_rev in searched:
                continue
            result = "{} {} {} {}".format(gene1, interaction_type, gene2, score)
            results.append(result)
            searched.append(query)
        return results

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
        try:
            check = self._db[query]
            return True
        except KeyError:
            return False

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
    >>> parsed_input = {"locus_name1" : ["gene1", "gene2", "gene3", "geneN"],
    >>>                 "locus_name2" : ["gene1", "gene2", "gene3", "geneN"]}
    """

    locus_genes = defaultdict(lambda: None)
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            data = line.strip("\n").split("\t")
            locus = data[1].split()[-1]
            genes = data[2:]
            locus_genes[locus] = genes

    return locus_genes
