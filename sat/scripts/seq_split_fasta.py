# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import make_output_dir, write_fasta, talk_to_me

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def seq_split_fasta(input_fasta: str, outfile_dir: str):
    """
    Creates single entry fasta files from a fasta file with multiple entries.
    """

    if input_fasta.endswith(".fasta"):
        write_fasta(input_fasta = input_fasta, outfile_dir = outfile_dir)

    else:
        talk_to_me("This is not a fasta file. Fasta file must be unzipped.")

    return 


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def seq_split_fasta_main(args):
    make_output_dir(args.outfile_dir, is_dir=True)
    seq_split_fasta(args.in_fasta, args.outfile_dir)
    return


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
