from collections import defaultdict
from datetime import datetime

def generate_id() -> str :
    """Generates a unique id based on time and date
    """
    unique_id = datetime.now().strftime("%M%d%y-%H%M%S")
    return unique_id


def save_as_sif(data : dict, out_dir="Simple_out") -> None:
    """ Converts data into sif format.
    Generates one SIF file per locus and also all loci in one SIF. The outputed files
    will contains a unique id that pertains to `Month Day Year - Hours Minutes Seconds'
    format. This allows the user to not only keep track when these files were generated
    but also to prevent overwriting.

    Arguments
    ---------
    data : dict
        Dictioanry containing the locus name and interaction
        scores
    out_dir : str (default="Simple_out")
        Name of the directory where all the generated SIF
        files will be stored

    Returns
    -------
    None
        Writes SIF files
    """
    # NOTE: add comment here
    unique_id = datetime.today().strftime("%m%d%y-%H%M%S")
    for locus_name, interactions in data.items():
        outname = "{}-{}.sif".format(locus_name, unique_id)
        with open(outname, "w") as sifile:
            sifile.write("gene1 interaction gene2 score\n")
            for interaction in interactions:
                sifile.write("{}\n".format(interaction))

def create_gmt_file(fa_data, networks, db, outname="Day3_Output"):
    """write out a gmt file

    Parameters
    ----------
    fa_data : dict
        FA locus data
    networks : list
        list of networks
    db : dict
        stringDB
    """
    locus_and_score = defaultdict(list)
    n_locus = len(fa_data.keys())
    for locus_idx in range(1, n_locus+1):
        fa_key = "locus {}".format(locus_idx)
        gene_arr = set(fa_data[fa_key])
        gene_score = defaultdict(list)
        sub_dict = defaultdict(list)
        for network in networks[:10]:
            selected_gene = network[locus_idx-1]
            locus_genes = list(gene_arr - set(selected_gene))
            for gene in locus_genes:
                network[locus_idx-1] = gene
                gene_score =  db.targeted_gene_scoring(gene, locus_genes)
                sub_dict[gene].append(gene_score)
        locus_and_score[fa_key]  = sub_dict

    scores = defaultdict(lambda: None)
    for locus, genes_scores in locus_and_score.items():
        sub_dict = defaultdict(lambda: None)
        for gene, score_arr in genes_scores.items():
            score = round(sum(score_arr)/len(networks),3)
            sub_dict[gene] = score
        scores[locus] = sub_dict

    # next write gmt file
    outfile_data = []
    for locus_id, genes in fa_data.items():
        scored_arr = scores[locus_id]
        line_out = []
        line_out.append(locus_id)
        for gene in genes:
            score = scored_arr[gene]
            if score == 0.0:
                score = "N/A"
                scored_gene = "{} {}".format(gene, score)
                line_out.append(scored_gene)
            else:
                scored_gene = "{} {}".format(gene, score)
                line_out.append(scored_gene)
        outfile_data.append(line_out)

    # converting into a string
    filename = "{}.txt".format(outname)
    with open(filename, "w") as outfile:
        for line_data in outfile_data:
            stdout = "\t".join(line_data)
            outfile.write("{}\n".format(stdout))