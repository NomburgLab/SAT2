import pytest
from sat.scripts.utils.misc import read_tsv


class Test_Read_Tsv():

    def test_header_from_file(self, tmp_path):
        """First line used as header when colnames not provided."""
        f = tmp_path / "data.tsv"
        f.write_text("name\tage\tcity\nAlice\t30\tNY\nBob\t25\tLA\n")

        rows = list(read_tsv(str(f)))
        assert len(rows) == 2
        assert rows[0] == {"name": "Alice", "age": "30", "city": "NY"}
        assert rows[1] == {"name": "Bob", "age": "25", "city": "LA"}

    def test_user_provided_colnames(self, tmp_path):
        """User-provided colnames, file has no header."""
        f = tmp_path / "data.tsv"
        f.write_text("Alice\t30\tNY\nBob\t25\tLA\n")

        rows = list(read_tsv(str(f), colnames="name,age,city"))
        assert len(rows) == 2
        assert rows[0] == {"name": "Alice", "age": "30", "city": "NY"}

    def test_user_provided_colnames_with_matching_header(self, tmp_path):
        """User provides colnames that match the first line (header is skipped)."""
        f = tmp_path / "data.tsv"
        f.write_text("name\tage\tcity\nAlice\t30\tNY\n")

        rows = list(read_tsv(str(f), colnames="name,age,city"))
        assert len(rows) == 1
        assert rows[0] == {"name": "Alice", "age": "30", "city": "NY"}

    def test_empty_file_raises(self, tmp_path):
        """Empty file raises ValueError."""
        f = tmp_path / "empty.tsv"
        f.write_text("")

        with pytest.raises(ValueError, match="empty"):
            list(read_tsv(str(f)))

    def test_blank_lines_skipped(self, tmp_path):
        """Blank lines in the middle and end are skipped."""
        f = tmp_path / "data.tsv"
        f.write_text("name\tage\n\nAlice\t30\n\nBob\t25\n\n")

        rows = list(read_tsv(str(f)))
        assert len(rows) == 2

    def test_only_blank_lines_raises(self, tmp_path):
        """File with only blank lines raises ValueError."""
        f = tmp_path / "blanks.tsv"
        f.write_text("\n\n\n")

        with pytest.raises(ValueError, match="empty"):
            list(read_tsv(str(f)))


