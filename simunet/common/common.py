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
