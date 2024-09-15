from .utils.misc import make_output_dir, talk_to_me
from .utils.ete3_taxonomy import Taxon

import pandas as pd


def validate_and_format_args(args):
    # Format args
    args.colnames = args.colnames.split(",")
    args.taxonomy_levels = args.taxonomy_levels.split(",")

    return args


def tab_add_taxonomy_main(args):
    args = validate_and_format_args(args)

    talk_to_me("Parsing infile")
    if args.colnames != [""]:
        df = pd.read_csv(args.infile, names=args.colnames, sep="\t")
    else:
        talk_to_me(
            "Colnames haven't been specified, so assuming the first row contains the "
            "colnames."
        )
        df = pd.read_csv(args.infile, header=0, sep="\t")

    # Check if 'taxid' is a column
    if "taxid" not in df.columns:
        msg = (
            f"One column must be named taxid. Currently, the colnames are {df.columns}."
            " (If this looks like a row of data, go back and add colnames.)"
        )
        raise ValueError(msg)

    lineage_levels = args.taxonomy_levels

    taxon_objects = dict()

    # Process each row
    talk_to_me("Adding taxonomy information")
    for index, row in df.iterrows():
        taxid = row["taxid"]

        if taxid in taxon_objects:
            taxon = taxon_objects[taxid]
        else:
            taxon = Taxon(taxid, lineage_levels)
            taxon_objects[taxid] = taxon

        canonical_lineage = taxon.canonical_lineage

        # Add lineage info to the DataFrame
        for level, name in zip(lineage_levels, canonical_lineage):
            df.at[index, level] = name

    # Write output
    talk_to_me("Writing output file")
    make_output_dir(args.outfile)
    df.to_csv(args.outfile, index=False, sep="\t")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
