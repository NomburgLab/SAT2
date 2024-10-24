from .utils.misc import talk_to_me, make_output_dir
from .utils.dali import parse_structure_key, read_alignment_block, segment_alignments, DALI_alignment


class DALI_alignment_motif_finder(DALI_alignment):

    def __init__(self, alignment, key):
        super().__init__(alignment, key)

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


def aln_dali_motif_finder_main(args):

    if args.key != "":
        key = parse_structure_key(args.key)
    else:
        key = ""
    alignments = read_alignment_block(args.aln_file)
    alignments = segment_alignments(alignments)

    out = ""
    total_aln_count = len(alignments)
    has_motif_count = 0
    for alignment in alignments:

        a = DALI_alignment_motif_finder()
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
