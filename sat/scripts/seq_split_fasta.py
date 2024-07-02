# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import make_output_dir, talk_to_me
import os
from Bio import SeqIO

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def seq_split_fasta(input_fasta: str, outfile_dir: str):
    """
    Creates single entry fasta files from a fasta file with multiple entries.
    """
    if input_fasta.endswith(".fasta"):
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            file_name = seq_record.id +'.fasta'
            file_path = os.path.join(outfile_dir, file_name) 
        
            with open(file_path, 'w') as file:
                file.write('>'+seq_record.description + '\n')
                file.write(str(seq_record.seq))

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
