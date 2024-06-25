# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import make_output_dir
import os
import gzip
from Bio import SeqIO


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def seq_split_fasta(input_fasta: str, outfile_dir: str):
    """
    Creates single entry fasta files from a fasta file with multiple entries.
    This function can handle .gzip files, but input_fasta needs to
    end with .gz
    """
   
    #Files are named after their accession. Maybe change naming?
    def write_fasta(input_fasta, outfile_dir):
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            file_name = seq_record.id +'.fasta'
            file_path = os.path.join(outfile_dir, file_name) 
            
            with open(file_path, 'w') as file:
                file.write('>'+seq_record.description + '\n')
                file.write(str(seq_record.seq))
        return


    if not input_fasta.endswith(".gz"):
        write_fasta(input_fasta = input_fasta, outfile_dir = outfile_dir)

    elif input_fasta.endswith(".gz"):
        with gzip.open(input_fasta, "rt") as handle:
            write_fasta(input_fasta = handle, outfile_dir = outfile_dir)
            

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
