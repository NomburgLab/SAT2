# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #

from .utils.Foldseek_Dataset import Foldseek_Dataset
from .utils.misc import read_fasta_to_memory, talk_to_me

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

def parse_aln_file(aln_dataset_object, aln_file, aln_fields):
    """
    Given an alignment file and the alignment fields/columns, create a dictionary of key:value pairs.
    The key is the target, and the value is a list of lists. Each list  contains the target start (tstart) and target end (tend) index,
    indicating the positions where the target aligns to the query. 

    Inputs:
    - aln_dataset_object: This is an instantiation of the Foldseek_Dataset class, which will store the alignment information
    in a tabular output file. 
    - aln_file: Alignment file from blast or mmseq. The file must contain the target, tstart,tend.
    - aln_fields: These are the column names for the alignment file. This is a comma delimited list of column names.

    Outputs:
    - aln_dict: A dictionary of targets and a list of lists, where each list is [tstart, tend, tlength]. 
      tlength is the length of the target that aligns to the query.
    """
    aln_dataset_object.parse_alignment(aln_file, aln_fields)

    aln_dict = {}

    for aln_key, aln_group in aln_dataset_object.alignment_groups.items():

        for aln_object in aln_group.alignments: 
    
            tstart = int(aln_object.tstart)
            tend = int(aln_object.tend)
            tlength = (tend - tstart)+1

            if aln_object.target not in aln_dict:
                aln_dict[aln_object.target] = []
                aln_dict[aln_object.target].append([tstart, tend, tlength])
            
            else:
                aln_dict[aln_object.target].append([tstart, tend, tlength])
    
    return aln_dict


def select_alignment(aln_list, short_aln_type = 1):
    """
    Based on user input, find the list of tstart, tend, tlength values
    for either the longest or shortest alignment length for the target.
    If there are multiple alignments that correspond to the longest or shortest alignment length,
    the first alignment that meets the requirement will be used.

    Inputs:
    - aln_list: A list of lists where each list contains tstart, tend, tlength (e.g. [tstart, tend, tlength]) for the target
    - short_aln_type: If short_aln_type is 1, find the shortest alignment length. 
                      If short_aln_type is 0, find the longest alignment length.
                      The default value is 1. 
    """
    short_aln_type = bool(short_aln_type)
    if short_aln_type == True:
        shortest_aln_len = min(aln_list, key=lambda x: x[2])
        return shortest_aln_len
    
    else:
        longest_aln_len = max(aln_list, key=lambda x: x[2])
        return longest_aln_len


def trim_sequence(chosen_aln, sequence):
    """
    Given a sequence and its corresponding list of tstart, tend, tlength values, trim the sequence using tstart and tend.
    tstart and tend must be 1-indexed.

    Inputs:
    - chosen_aln: A list of tstart, tend, tlength values that corresponds to either the longest or shortest alignment.
    - sequence: protein sequence

    Outputs:
    - trimmed_sequence: A trimmed sequence from position tstart to position tend.

    """
    if len(chosen_aln)!=3:
        raise ValueError(f"The alignment list for the target should be [tstart, tend, tlength]. \
                         {chosen_aln} doesn't match the expected format.")
    
    tstart = chosen_aln[0] - 1
    tend = chosen_aln[1] - 1
    tlength = chosen_aln[2]
    trimmed_sequence = sequence[tstart:tend+1]

    if len(trimmed_sequence) != tlength:
        raise ValueError(f"The length of trimmed target sequence does not match its expected length.")
    
    return trimmed_sequence
    
def aln_trim_target_fasta_main(args):
    """
    Generate a fasta file of target accessions and their trimmed sequences. 
    
    Inputs:
    - alignment_file: Alignment file from blast or mmseq. The file must contain the target, tstart,tend.
    - alignment_fields: These are the column names for the alignment file. This is a comma delimited list of column names.
    - fasta_file: Fasta file of the targets and their sequences
    - output_file: name for output file
    - short_aln_type: Indicate 1 for shortest alignment and 0 for longest alignment for the targets

    Outputs:
    - output_file: a fasta file with the target accessions and their corresponding trimmed sequences
    """

    #check alignment fields for tstart and tend 
    if 'tstart' and 'tend' not in args.alignment_fields:
        raise ValueError(f"tstart and tend alignment fields are not found in alignment_fields parameter.")

    aln_dataset_object = Foldseek_Dataset()
    aln_dict = parse_aln_file(aln_dataset_object, args.alignment_file, args.alignment_fields)
    fasta_dict = read_fasta_to_memory(args.fasta_file)

    #for each target, output its trimmed sequence to the output file
    output_file = open(args.output_file, "a")

    for target, aln_list in aln_dict.items():
        chosen_aln = select_alignment(aln_list, args.short_aln_type)

        if target not in fasta_dict:
            talk_to_me(f"Target accession '{target}' not found in the fasta file.")

        else:
            target_sequence = fasta_dict[target]
            trimmed_target_sequence = trim_sequence(chosen_aln, target_sequence)

            output_file.write('>' + target + "\n")
            output_file.write(trimmed_target_sequence + "\n")
    
    output_file.close()

if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
