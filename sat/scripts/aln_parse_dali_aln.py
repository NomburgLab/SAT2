import re

from .utils.misc import talk_to_me, make_output_dir
from .utils.dali import parse_structure_key


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

    def parse_alignment(self, alignment, key):
        """
        This takes in a single alignment block, and extracts the alignment information
        (including DSSPs, qseq, tseq, etc.). Also calls parse_alignment_header to parse
        the header, using the dali key.

        The key variable can be empty - then, just doesn't fill out the query and
        target names (although query_id and target_id will be there still)
        """

        # Parse the header, then toss it
        self.parse_alignment_header(alignment[0], key)
        alignment = alignment[1:]

        aln_qDSSP = ""
        aln_qseq = ""
        aln_tseq = ""
        aln_tDSSP = ""

        for i, line in enumerate(alignment):

            if line.startswith("DSSP"):
                line = line[6:]
                if alignment[i - 1].startswith("Sbjct"):
                    aln_tDSSP += line
                elif alignment[i + 1].startswith("Query"):
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
            print(alignment)
            raise ValueError(msg)

        self.aln_qDSSP = aln_qDSSP
        self.aln_qseq = aln_qseq
        self.aln_tseq = aln_tseq
        self.aln_tDSSP = aln_tDSSP

    def parse_alignment_header(self, header, key):
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

    def get_alignment_position(self, query_position):
        """
        This function takes in the residue position of a query, and returns where along
        the alignment that residue is located.

        Notes on indexing:
        - query_position is 1-indexed. E.g. 4 indicates the 4th residue.
        - The output value is 0-indexed to allow easy subsetting later on.
        """
        count = 0
        for i, field in enumerate(self.aln_qseq):
            if field != "-":
                count += 1
            if count == query_position:
                return i

        # if we get here, there is an error
        msg = (
            "Erorr in get_alignment_position()! The input query_position is "
            f"{query_position}, which is a {type(query_position)}. The self.aln_qseq is {self.aln_qseq}"
        )
        raise ValueError(msg)

    def alignment_has_residues(self, aln_pos, residues):
        """
        This function takes in a 1-indexed residue position and a more complicated
        'residues' string that indicates which residues in the target are desired in that
        alignment position..

        aln_pos: The 0-indexed alignment position at which the original query position
            falls.
        residues: This is a comma-delimited string of residues. For each item in this list,
            however, all possible residues are separated by forward slashes. If a
            residue is listed as 'X', this means that any residue (or alignment gap) is
            fine.

        Example:
        - position: 5
        - residues: b/z,c/t,d/q

        Result:
        - This will look up where in the alignment the 1-indexed 5th residue of the query
        is.
        - It will then determine that that alginment position is a B or Z in the target,
        followed by a C or T, followed by a D or Q.
        """
        try:
            aln_pos = int(aln_pos)
        except ValueError:
            msg = f"The 'aln_pos' input here, {aln_pos}, isn't an integer!"
            raise ValueError(msg)

        residues = residues.upper().split(",")
        motif_len = len(residues)

        for i in range(motif_len):

            residue = residues[i]
            allowed_residues = residue.split("/")

            current_pos = aln_pos + i

            if allowed_residues == ["X"]:
                continue

            if self.aln_tseq[current_pos].upper() in allowed_residues:
                continue
            else:
                return False

        # if we made it, it is a success
        return True

    def has_motif(self, motif):
        """
        The motif is in format of pos_residues_flexibility, as per the follows:
        - pos: 1-indexed integer value of a specific integer in the query. This code
            will call get_alignment_position() to determine where in the alignment
            this residue falls.
        - residues: a comma-delimited list of residues. The length of that list is the
            length of the motif - e.g. H,H means you are looking for two sequential
            histidines in the target at the alignment position. However, multiple residue
            possibilities can be indicated for each position. E.g. H,T/S indicates we are
            looking for an H followed by a T or S. The alignment position indicates the
            start of this motif. You could also do something like H,X,H which means
            histidine followed by any residue (or an alignment gap) followed by another
            histidine.
        - flexibility: this is where the motif can start relative to the indicated pos. A
            value of 0 means that the motif should start exactly at the alignment position
            derived by pos. A value of 1 means that the motif can start 1 residue before
            or 1 residue after the derived alignment position.

        Some valid examples of the motif input:
        - 2_b/z,c/t,x_1
        - 2_B/Z,C/T,X_1
        - 15_H_0
        - etc...
        """
        position, residues, flexibility = motif.split("_")
        position = int(position)
        flexibility = int(flexibility)
        aln_pos = self.get_alignment_position(position)
        for i in range(aln_pos - flexibility, aln_pos + flexibility + 1):
            if self.alignment_has_residues(i, residues):
                return True
        return False


def aln_parse_dali_aln_main(args):

    key = parse_structure_key(args.key)
    alignments = read_alignment_block(args.aln_file)
    alignments = segment_alignments(alignments)

    out = ""
    total_aln_count = len(alignments)
    has_motif_count = 0
    for alignment in alignments:

        a = DALI_alignment()
        a.parse_alignment(alignment, key)

        has_all_motifs = True
        for motif in args.motif_list.split("+"):
            if not a.has_motif(motif):
                has_all_motifs = False
                break

        if has_all_motifs:
            has_motif_count += 1
            out += f"{a.target}\t{a.target_id}\n"

    talk_to_me(f"{has_motif_count} of {total_aln_count} targets have all motifs.")

    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
