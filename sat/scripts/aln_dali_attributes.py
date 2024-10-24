import math
import pandas as pd
#from .aln_parse_dali_matrix import parse_key
from .utils.misc import talk_to_me
from .utils.dali import DALI_alignment, read_alignment_block, segment_alignments, parse_structure_key


class DALI_alignment_attributes(DALI_alignment):

    def __init__(self, alignment, key = ""):
        super().__init__(alignment, key)
        self.parse_alignment()

    def find_aln_position(self):
        """
        Using the query and target alignment sequences from a DALI alignment block, determine the first and last indices of their alignment.
        Residues only align if they are upper case in both the query and target alignment sequences, at the same index. 
        Indices containing hyphens or lower case residues do not incidate alignment.
        The alignments are zero-indexed.

        Inputs:
        - aln_qseq: This is the query alignment sequence from a DALI alignment object, containing hyphens and amino acids.
        - aln_tseq: This is the target alignment sequence from a DALI alignment object, containing hyphens and amino acids.
        
        Outputs:
        - aln_start: This is the first position of the alignment block where the query and target align.
        - aln_end: This is the last position of the alignment block where the query and target align.
        """
        #check if input attributes exist
        if not hasattr(self, 'aln_qseq'):
            raise ValueError("Alignment query sequence attribute not found.")
        if not hasattr(self, 'aln_tseq'):
            raise ValueError("Alignment target sequence attribute not found.")

        aln_pos_list = []
        aln_index = -1

        #loop through the query and target alignment sequences and find structurally equivalent residues, in which both amino acids are upper case
        for query_aa, target_aa in zip(self.aln_qseq, self.aln_tseq):
            aln_index  +=1

            if query_aa!='-' and target_aa!='-':

                #Alignment occurs when the amino acids of the target and query are both upper case 
                if query_aa.isupper() and target_aa.isupper():
                    aln_pos_list.append(aln_index)
                
                #check if the amino acid of the target or query is upper case and the other is lower case. This shouldn't happen.
                if query_aa.isupper() ^ target_aa.isupper():
                    raise ValueError("One amino acid in the query or target is uppercase/aligned, while the other is not!")

        if len(aln_pos_list) == 0:
            raise ValueError("No valid alignment positions found.")
        
        #set aln_start and aln_end as attributes
        aln_start = aln_pos_list[0]
        aln_end = aln_pos_list[-1]

        self.aln_start = aln_start
        self.aln_end = aln_end

    def find_seq_position(self, aln_position, aln):
        """
        An alignment, whether for the query or target, consists of both hyphens and amino acids.
        The hyphens and lower case amino acids indicates where the query does not align with the target (or vice versa). 
        Given an index position in the alignment, determine the corresponding position in the protein sequence.
        The protein sequence only contains residues/amino acids.

        Inputs:
        - aln_position: Index of the alignment position. The alignment is zero-indexed.
        - aln: The alignment string that includes both hyphens and amino acids.

        Outputs:
        - seq_position: Index of the amino acid in the protein sequence that corresponds to the specified alignment position.
                        The protein sequence is 1-indexed.
        """

        if aln_position >= len(aln):
            raise ValueError("Alignment position is greater than the alignment length.")
        
        index_count = 0

        for aln_index, aa in enumerate(aln):
            if aa != "-":
                index_count+=1
                
                if aln_index == aln_position:
                    return index_count
        
        raise ValueError("Sequence position not found.")

    
    def calculate_pident(self): 
        """
        pident refers to %id in DaliLite.v5 for structural alignments.
        %id is the percentage of identical amino acids out of structurally equivalent residues.
        Structurally equivalent residues are upper case.

        Inputs:
        - aln_qseq: This is the query alignment sequence from a DALI alignment object, containing hyphens and amino acids.
        - aln_tseq: This is the target alignment sequence from a DALI alignment object, containing hyphens and amino acids.

        Outputs:
        - pident: Pident is percentage of identical amino acids out of structurally equivalent residues. 
                  Pident is rounded to the nearest whole integer. 
        """
        #check if input attributes exist
        if not hasattr(self, 'aln_qseq'):
            raise ValueError("Alignment query sequence attribute not found.")
        if not hasattr(self, 'aln_tseq'):
            raise ValueError("Alignment target sequence attribute not found.")
        
        #create a custom function to always round up when the decimal has 0.5. 
        def custom_round(x):
            frac = x - math.floor(x)
            if frac < 0.5: 
                return math.floor(x)
            return math.ceil(x)
        
        #count the number of structurally equivalent residues and total number of aligned residues to calculate pident
        identical_aa_count = 0
        aligned_residues = 0  
        
        for query_aa, target_aa in zip(self.aln_qseq, self.aln_tseq):

            if query_aa != '-' and target_aa != '-':

                #check for alignment. Amino acids that are upper case are aligned
                if query_aa.isupper() and target_aa.isupper():    
                    aligned_residues += 1

                    #check for identical residues b/w query and targety
                    if query_aa == target_aa:
                        identical_aa_count += 1  
        
        #calculate pident. We will need to round pident to a whole number.
        if aligned_residues == 0:
            raise ValueError("Alignment query and target have no aligned residues...pident will not be calculated.")
        
        pident = (identical_aa_count / aligned_residues) * 100
        pident = custom_round(pident) 

        return pident


    def get_attribute_dict(self):
        """
        Store the attributes of this query/target alignment in a dictionary.

        Outputs:
        - aln_dict: A dictionary of attributes for a query/target alignment.
                    Attributes include query_id, target_id, qstart, qend, tstart, etc.
        """
        #check if input attributes exist
        check_attr_list = ["aln_start", "aln_qseq", "aln_end", "aln_qseq", "aln_tseq"]
        for attr in check_attr_list:
            if not hasattr(self, attr):
                raise ValueError(f"The {attr} attribute is not found.")
        
        #determine start and end indices of query protein sequence
        self.qstart = self.find_seq_position(self.aln_start, self.aln_qseq)
        self.qend = self.find_seq_position(self.aln_end, self.aln_qseq)
        
        #determine start and end indices of target protein sequence
        self.tstart = self.find_seq_position(self.aln_start, self.aln_tseq)
        self.tend = self.find_seq_position(self.aln_end, self.aln_tseq)
        
        #calculate pident
        self.pident  = self.calculate_pident()

        #create a dictionary with all the attributes
        aln_dict = {}
        aln_dict['aln_number'] = self.aln_number
        aln_dict['query'] = self.query
        aln_dict['target'] = self.target
        aln_dict['qstart'] = self.qstart
        aln_dict['qend'] = self.qend
        aln_dict['tstart'] = self.tstart
        aln_dict['tend'] = self.tend
        aln_dict['qlen'] = self.seq_qlen
        aln_dict['tlen'] = self.seq_tlen
        aln_dict['z'] = self.z
        aln_dict['aln_number'] = self.aln_number 
        aln_dict['pident'] = self.pident 
        
        return aln_dict

def aln_dali_alignment_attributes_main(args):
    """
    Generates a csv file of aligned targets and queries and their attributes.
    Attributes such as qstart, qend, tstart, tend, etc are retrieved from the alignment object. 
    
    Inputs:
    - aln_file: The alignment file is a text file output by DALI containing pairwise alignments, summary block, etc.
    - key: This is a text file of the DALI database key. 
            The key contains the structure file name and its associated 4 digit DALI identifier. 
    - outfile: Path to the csv outfile.

    Outputs:
    - output_file: a csv file containing alignment pairs and their attributes.
    """
        
    contents = read_alignment_block(args.aln_file) 
    segmented_aln = segment_alignments(contents)

    aln_dict_list  = []

    if args.key != "":
        key = parse_structure_key(args.key)

    #for each alignment block, determine its attributes. Store the attributes in a dictionary. 
    for aln in segmented_aln:
        aln_obj = DALI_alignment_attributes(aln, key)
        aln_obj.find_aln_position()
        aln_obj.info = aln_obj.get_attribute_dict()

        #add the dictionary of attributes to a list
        aln_dict_list.append(aln_obj.info)
    
    #convert the list of attribute dictionaries and output it as a csv file
    df = pd.DataFrame(aln_dict_list)
    df.to_csv(args.outfile, index=False)