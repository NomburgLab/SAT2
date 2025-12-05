import pytest
import tempfile
import os

from sat.scripts.aln_merge_clusters import (
    remove_file_suffix,
    parse_cluster_file,
    merge_cluster_files,
)


# ------------------------------------------------------------------------------------ #
# Test data paths
# ------------------------------------------------------------------------------------ #
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data", "foldseek_related")
FILE1_SIMPLE = os.path.join(TEST_DATA_DIR, "merge_clusters_file1_simple.tsv")
FILE2_NESTED = os.path.join(TEST_DATA_DIR, "merge_clusters_file2_nested.tsv")
FILE2_SIMPLE = os.path.join(TEST_DATA_DIR, "merge_clusters_file2_simple.tsv")


# ------------------------------------------------------------------------------------ #
# Tests for remove_file_suffix
# ------------------------------------------------------------------------------------ #
def test_remove_file_suffix__pdb():
    assert remove_file_suffix("protein.pdb") == "protein"


def test_remove_file_suffix__fasta():
    assert remove_file_suffix("sequence.fasta") == "sequence"


def test_remove_file_suffix__fa():
    assert remove_file_suffix("sequence.fa") == "sequence"


def test_remove_file_suffix__cif():
    assert remove_file_suffix("structure.cif") == "structure"


def test_remove_file_suffix__pdb_gz():
    assert remove_file_suffix("protein.pdb.gz") == "protein"


def test_remove_file_suffix__no_suffix():
    assert remove_file_suffix("protein") == "protein"


def test_remove_file_suffix__case_insensitive():
    assert remove_file_suffix("protein.PDB") == "protein"
    assert remove_file_suffix("sequence.FASTA") == "sequence"


def test_remove_file_suffix__faa():
    assert remove_file_suffix("sequence.faa") == "sequence"


# ------------------------------------------------------------------------------------ #
# Tests for parse_cluster_file
# ------------------------------------------------------------------------------------ #
def test_parse_cluster_file__with_header():
    colnames, rows = parse_cluster_file(FILE1_SIMPLE)
    assert colnames == ["lol_rep", "struc_rep"]
    assert len(rows) == 3
    # Check that suffixes are removed
    assert rows[0]["lol_rep"] == "A"
    assert rows[0]["struc_rep"] == "A"


def test_parse_cluster_file__nested():
    colnames, rows = parse_cluster_file(FILE2_NESTED)
    assert colnames == ["struc_rep", "seq_rep", "member"]
    assert len(rows) == 5
    # Check specific row
    assert rows[0]["struc_rep"] == "A"
    assert rows[0]["seq_rep"] == "B"
    assert rows[0]["member"] == "B"


def test_parse_cluster_file__with_provided_colnames():
    """Test providing column names explicitly."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("X.pdb\tY.pdb\n")
        f.write("A.pdb\tB.pdb\n")
        temp_path = f.name

    try:
        colnames, rows = parse_cluster_file(temp_path, "col1,col2")
        assert colnames == ["col1", "col2"]
        assert len(rows) == 2
        assert rows[0]["col1"] == "X"
        assert rows[0]["col2"] == "Y"
    finally:
        os.unlink(temp_path)


def test_parse_cluster_file__colnames_match_header():
    """Test when provided colnames match the first line (should skip header)."""
    colnames, rows = parse_cluster_file(FILE1_SIMPLE, "lol_rep,struc_rep")
    assert colnames == ["lol_rep", "struc_rep"]
    assert len(rows) == 3  # Should not include header as data


def test_parse_cluster_file__wrong_column_count():
    """Test error when line has wrong number of columns."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("col1\tcol2\n")
        f.write("A\tB\tC\n")  # Too many columns
        temp_path = f.name

    try:
        with pytest.raises(ValueError, match="columns but expected"):
            parse_cluster_file(temp_path)
    finally:
        os.unlink(temp_path)


# ------------------------------------------------------------------------------------ #
# Tests for merge_cluster_files
# ------------------------------------------------------------------------------------ #
def test_merge_cluster_files__simple_to_nested():
    """Test merging a simple file1 with a nested file2."""
    file1_colnames, file1_rows = parse_cluster_file(FILE1_SIMPLE)
    file2_colnames, file2_rows = parse_cluster_file(FILE2_NESTED)

    out_colnames, out_rows, unmatched_count = merge_cluster_files(
        file1_colnames, file1_rows, file2_colnames, file2_rows
    )

    # Check output columns
    assert out_colnames == ["lol_rep", "struc_rep", "seq_rep", "member"]

    # Check output row count
    assert len(out_rows) == 5

    # Check no unmatched entries
    assert unmatched_count == 0

    # Check specific mappings
    # Row with member B should have lol_rep=A (since struc_rep A maps to lol_rep A)
    b_rows = [r for r in out_rows if r["member"] == "B"]
    assert len(b_rows) == 1
    assert b_rows[0]["lol_rep"] == "A"
    assert b_rows[0]["struc_rep"] == "A"
    assert b_rows[0]["seq_rep"] == "B"

    # Row with member D should have lol_rep=A (since struc_rep D maps to lol_rep A)
    d_rows = [r for r in out_rows if r["member"] == "D"]
    assert len(d_rows) == 1
    assert d_rows[0]["lol_rep"] == "A"
    assert d_rows[0]["struc_rep"] == "D"

    # Row with member F should have lol_rep=E
    f_rows = [r for r in out_rows if r["member"] == "F"]
    assert len(f_rows) == 1
    assert f_rows[0]["lol_rep"] == "E"


def test_merge_cluster_files__simple_to_simple():
    """Test merging two simple 2-column files."""
    file1_colnames, file1_rows = parse_cluster_file(FILE1_SIMPLE)
    file2_colnames, file2_rows = parse_cluster_file(FILE2_SIMPLE)

    out_colnames, out_rows, unmatched_count = merge_cluster_files(
        file1_colnames, file1_rows, file2_colnames, file2_rows
    )

    # Check output columns
    assert out_colnames == ["lol_rep", "struc_rep", "member"]

    # Check output row count
    assert len(out_rows) == 5

    # Check no unmatched entries
    assert unmatched_count == 0


def test_merge_cluster_files__missing_in_file1():
    """Test when file2's first column value not found in file1 - uses its own value."""
    file1_colnames = ["lol_rep", "struc_rep"]
    file1_rows = [{"lol_rep": "A", "struc_rep": "A"}]

    file2_colnames = ["struc_rep", "member"]
    file2_rows = [
        {"struc_rep": "A", "member": "A"},
        {"struc_rep": "Z", "member": "Z"},  # Z not in file1
    ]

    out_colnames, out_rows, unmatched_count = merge_cluster_files(
        file1_colnames, file1_rows, file2_colnames, file2_rows
    )

    # Should have 1 unmatched entry
    assert unmatched_count == 1

    # Check the matched row
    a_rows = [r for r in out_rows if r["member"] == "A"]
    assert len(a_rows) == 1
    assert a_rows[0]["lol_rep"] == "A"

    # Check the unmatched row uses its own value as higher-level rep
    z_rows = [r for r in out_rows if r["member"] == "Z"]
    assert len(z_rows) == 1
    assert z_rows[0]["lol_rep"] == "Z"  # Uses Z as its own higher-level rep
    assert z_rows[0]["struc_rep"] == "Z"


def test_merge_cluster_files__missing_in_file2():
    """Test when file1's second column value not found in file2 - should be ignored."""
    file1_colnames = ["lol_rep", "struc_rep"]
    file1_rows = [
        {"lol_rep": "A", "struc_rep": "A"},
        {"lol_rep": "A", "struc_rep": "Z"},  # Z not in file2 - this is fine
    ]

    file2_colnames = ["struc_rep", "member"]
    file2_rows = [{"struc_rep": "A", "member": "A"}]

    out_colnames, out_rows, unmatched_count = merge_cluster_files(
        file1_colnames, file1_rows, file2_colnames, file2_rows
    )

    # Should work fine - unused file1 entries are just ignored
    assert len(out_rows) == 1
    assert unmatched_count == 0
    assert out_rows[0]["lol_rep"] == "A"


def test_merge_cluster_files__multiple_unmatched():
    """Test multiple unmatched entries from file2."""
    file1_colnames = ["lol_rep", "struc_rep"]
    file1_rows = [{"lol_rep": "A", "struc_rep": "A"}]

    file2_colnames = ["struc_rep", "member"]
    file2_rows = [
        {"struc_rep": "A", "member": "A"},
        {"struc_rep": "X", "member": "X"},  # Not in file1
        {"struc_rep": "Y", "member": "Y"},  # Not in file1
        {"struc_rep": "Z", "member": "Z"},  # Not in file1
    ]

    out_colnames, out_rows, unmatched_count = merge_cluster_files(
        file1_colnames, file1_rows, file2_colnames, file2_rows
    )

    # Should have 3 unmatched entries
    assert unmatched_count == 3
    assert len(out_rows) == 4

    # Check all unmatched rows use their own struc_rep value as lol_rep
    for row in out_rows:
        if row["member"] != "A":
            # lol_rep should equal struc_rep (uses its own value as higher-level rep)
            assert row["lol_rep"] == row["struc_rep"]


def test_merge_cluster_files__inconsistent_file1():
    """Test error when file1 has inconsistent mappings."""
    file1_colnames = ["lol_rep", "struc_rep"]
    file1_rows = [
        {"lol_rep": "A", "struc_rep": "X"},
        {"lol_rep": "B", "struc_rep": "X"},  # X maps to both A and B
    ]

    file2_colnames = ["struc_rep", "member"]
    file2_rows = [{"struc_rep": "X", "member": "X"}]

    with pytest.raises(ValueError, match="Inconsistent mapping"):
        merge_cluster_files(file1_colnames, file1_rows, file2_colnames, file2_rows)
