import re

def parse_structure_key(structure_key_file, delim=",,", existing_dict=dict()):
    """
    Takes in a path to a file with the first column the structure, and second column
    the 4-digit identifier, and returns a dictionary of format
    identifier:structure.

    Add an existing_dict if you just want to update that dictionary - this lets you
    have multiple structure_key files.
    """

    key_to_structure = existing_dict.copy()
    with open(structure_key_file) as infile:
        for line in infile:
            line = line.rstrip("\n")
            structure, key = line.split(delim)

            if len(key) != 4:
                msg = "The key is expected to be a four-digit identifier! This key, "
                msg += f" {key}, is a length of {len(key)}!"
                raise ValueError(msg)

            if key in key_to_structure:
                msg = f"Have obseved a key, {key}, that is already present in"
                msg += " key_to_structure! This means it may be present multiple times!"
                raise ValueError(msg)
            
            # unfortuantely my structure keys end in .pdb, which is not needed.
            structure = structure.rstrip(".pdb")

            key_to_structure[key] = structure

    return key_to_structure

def read_alignment_block(aln_file):
    """
    This function reads through a DALI output file and parses the lines that belong
    to the "Pairwise alignments" section of the file. It ignores the rest.

    Each line is then saved as an item in a list, which is the output of this function.
    """

    IN_ALN = False
    contents = []

    with open(aln_file) as infile:
        for line in infile:

            # Start of the alignment block
            if line.startswith("# Pairwise alignments"):
                IN_ALN = True
                continue

            # End of the alignment block. If we are currently in the alignment block but
            # now it is ending, can break the loop
            elif line.startswith("#"):
                if IN_ALN:
                    break

            if IN_ALN:
                line = line.rstrip("\n")
                if line == "":
                    continue
                contents.append(line)

    if contents == []:
        msg = (
            "Never found an alignment block in the input file! The alignment block"
            " should start with '# Pairwise alignments'"
        )
        raise ValueError(msg)

    return contents


def segment_alignments(alignments):
    """
    The input, alignments, is a list with each item being a line from the alignment
    chunk of a dali file - provided by read_alignment_block().

    In the alignment block, you have a bunch of different query-target pairs all
    concatenated one after another. This function splits the input into a list of lists,
    where each sublist contains the lines for a single query-target pair.
    """
    segmented_alignments = []
    current_segment = []

    for line in alignments:

        if line.startswith("No"):

            # Write the previous segment. We're now starting a new one
            if current_segment != []:
                segmented_alignments.append(current_segment)
                current_segment = []
        current_segment.append(line)

    # The last segment needs to be appended
    segmented_alignments.append(current_segment)
    return segmented_alignments


class DALI_alignment:

    def __init__(self, alignment, key = "") -> None:
        self.alignment = alignment
        self.key = key

    def parse_alignment(self):
        """
        This takes in a single alignment block, and extracts the alignment information
        (including DSSPs, qseq, tseq, etc.). Also calls parse_alignment_header to parse
        the header, using the dali key.

        The key variable can be empty - then, just doesn't fill out the query and
        target names (although query_id and target_id will be there still)
        """

        # Parse the header, then toss it
        self.parse_alignment_header(self.alignment[0], self.key)
        self.alignment = self.alignment[1:]

        aln_qDSSP = ""
        aln_qseq = ""
        aln_tseq = ""
        aln_tDSSP = ""

        for i, line in enumerate(self.alignment):

            if line.startswith("DSSP"):
                line = line[6:]
                if self.alignment[i - 1].startswith("Sbjct"):
                    aln_tDSSP += line
                elif self.alignment[i + 1].startswith("Query"):
                    aln_qDSSP += line

            elif line.startswith("Query"):
                line = line[6:].replace(" ", "")
                line = re.sub(r"\d+", "", line)  # removes trailing number
                aln_qseq += line

            elif line.startswith("Sbjct"):
                line = line[6:].replace(" ", "")
                line = re.sub(r"\d+", "", line)  # removes trailing number
                aln_tseq += line

        if not (len(aln_qDSSP) == len(aln_qseq) == len(aln_tseq) == len(aln_tseq)):
            msg = "When parsing the alignments for this chunk, got different lengths of one of"
            msg += " the lines."
            print(self.alignment)
            raise ValueError(msg)

        self.aln_qDSSP = aln_qDSSP
        self.aln_qseq = aln_qseq
        self.aln_tseq = aln_tseq
        self.aln_tDSSP = aln_tDSSP

        self.seq_qlen = self.calc_seq_len(self.aln_qseq)
        self.seq_tlen = self.calc_seq_len(self.aln_tseq)
    

    def parse_alignment_header(self, header, key ):
        """
        Parse the header part of an alignment chunk. Also takes in the key, to convert
        query and subject.

        Input:
        - 'No 3: Query=8IXZA Sbjct=cwtuA Z-score=7.8'

        Output:
        - Fills out the .query, .query_id, .target, .target_id, .z, and aln_number slots
        of self.
        """
        if not header.startswith("No "):
            msg = (
                "This function is supposed to parse the header part of an alignment chunk"
                " - e.g. 'No 3: Query=8IXZA Sbjct=cwtuA Z-score=7.8'"
            )
            raise ValueError(msg)

        # Parse the header
        parts = header.split()
        aln_number = parts[1].strip(":")
        query_id = parts[2].split("=")[1][:-1]
        target_id = parts[3].split("=")[1][:-1]
        z = float(parts[4].split("=")[1])

        # Convert query and subject if key provided
        query = "x"
        target = "x"
        if key != "":
            query = key.get(query_id, "x")[:-4]
            target = key.get(target_id, "x")[:-4]

        # Load self
        self.query = query
        self.query_id = query_id
        self.target = target
        self.target_id = target_id
        self.z = z
        self.aln_number = aln_number

    def __str__(self):
        ID = f"{self.query_id}_{self.target_id}"
        return ID

    def __repr__(self):
        ID = f"{self.query_id}_{self.target_id}"
        return ID

    def calc_seq_len(self, seq):
        seq_cleaned = seq.replace("-", "")
        if len(seq_cleaned) == 0:
            raise ValueError("Sequence does not contain amino acids, just hyphens.")
        return len(seq_cleaned) 