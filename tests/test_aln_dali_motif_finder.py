from sat.scripts.utils.dali import DALI_alignment
from sat.scripts.aln_dali_motif_finder import DALI_alignment_motif_finder


class Test_get_alignment_position:
    def test_1(self):
        a = DALI_alignment_motif_finder("", "")

        a.aln_qseq = "abcde--fghi--JK"
        pos = 5

        result = a.get_alignment_position(pos)

        assert result == 4
        assert a.aln_qseq[result] == "e"

    def test_2(self):
        a = DALI_alignment_motif_finder("", "")

        a.aln_qseq = "abcde--fghi--JK"
        pos = 6

        result = a.get_alignment_position(pos)

        assert result == 7
        assert a.aln_qseq[result] == "f"

    def test_3(self):
        a = DALI_alignment_motif_finder("","")

        a.aln_qseq = "abcde--fghi--JK"
        pos = 11

        result = a.get_alignment_position(pos)

        assert result == 14
        assert a.aln_qseq[result] == "K"


class Test_alignment_has_residues:
    def test_1(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "abcde--fghi--JK"

        motif = "2_b/z,c/t,d/q_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert self.alignment_has_residues(aln_position, residues)

    def test_1_uppercase_motif(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "abcde--fghi--JK"

        motif = "2_B/Z,C/T,D/Q_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert self.alignment_has_residues(aln_position, residues)

    def test_2(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "a-cde--fghi--JK"

        motif = "2_b/z,c/t,d/q_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert not self.alignment_has_residues(aln_position, residues)

    def test_3(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---cedjdafn-fds"

        motif = "7_a/b,f_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert self.alignment_has_residues(aln_position, residues)

    def test_3_with_X(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---cedjdafn-fds"

        motif = "7_a/b,f,x_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert self.alignment_has_residues(aln_position, residues)

    def test_3_with_X_FALSE(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---cedjdafn-fds"

        motif = "7_x,t,x_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert not self.alignment_has_residues(aln_position, residues)

    def test_single(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "--heqqew--sdfs-"

        motif = "5_q_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert self.alignment_has_residues(aln_position, residues)

    def test_single_false(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "--hexqew--sdfs-"

        motif = "5_q_0"
        position, residues, flexability = motif.split("_")
        position = int(position)
        aln_position = self.get_alignment_position(position)

        assert not self.alignment_has_residues(aln_position, residues)


class Test_has_motif:
    def test_has_motif_no_flex_true(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---cedjdafn-fds"

        motif = "7_a/b,f_0"

        assert self.has_motif(motif)

    def test_has_motif_no_flex_false(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---cedjdafn-fds"

        motif = "7_q,f_0"

        assert not self.has_motif(motif)

    def test_has_motif_flex_true(self):
        self = DALI_alignment_motif_finder("","")
        self.aln_qseq = "abcde--fghi--JK"
        self.aln_tseq = "---de-fghfn-fds"

        motif = "6_f,g,h_1"

        assert self.has_motif(motif)
