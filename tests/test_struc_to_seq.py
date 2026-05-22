import os

import pytest

from sat.scripts.struc_to_seq import struc_to_seq_main, _extract_seq
from sat.scripts.utils.structure import pdb_to_structure_object, struc_to_seq


PDB_FILE = "tests/test_data/structure_related/discontinuous_structure.pdb"
EXPECTED_SEQ = struc_to_seq(pdb_to_structure_object(PDB_FILE))


# ------------------------------------------------------------------------------------ #
# Unit tests for _extract_seq
# ------------------------------------------------------------------------------------ #
def test_extract_seq_ok():
    seq_record = _extract_seq(PDB_FILE)
    assert seq_record.startswith(">discontinuous_structure\n")
    seq = seq_record.strip().split("\n")[1]
    assert seq == EXPECTED_SEQ


def test_extract_seq_bad_file(tmp_path):
    bad_pdb = tmp_path / "bad.pdb"
    bad_pdb.write_text("not a pdb file\n")
    with pytest.raises(ValueError, match=str(bad_pdb)):
        _extract_seq(str(bad_pdb))


# ------------------------------------------------------------------------------------ #
# Single-file mode tests
# ------------------------------------------------------------------------------------ #
def test_single_file_print(capsys):
    class args:
        structure_file = PDB_FILE
        out_file = ""
        header = ""

    struc_to_seq_main(args)
    captured = capsys.readouterr()
    assert EXPECTED_SEQ in captured.out


def test_single_file_write(tmp_path):
    out = str(tmp_path / "out.fasta")

    class args:
        structure_file = PDB_FILE
        out_file = out
        header = "test_header"

    struc_to_seq_main(args)

    with open(out) as f:
        content = f.read()
    assert ">test_header\n" in content
    assert EXPECTED_SEQ in content


def test_single_file_append(tmp_path):
    out = str(tmp_path / "out.fasta")

    for hdr in ["seq1", "seq2"]:
        class args:
            structure_file = PDB_FILE
            out_file = out
            header = hdr
        struc_to_seq_main(args)

    with open(out) as f:
        content = f.read()
    assert content.count(">") == 2
    assert ">seq1\n" in content
    assert ">seq2\n" in content


# ------------------------------------------------------------------------------------ #
# Directory mode tests
# ------------------------------------------------------------------------------------ #
def _make_test_dir(tmp_path):
    """Copy test PDBs into a temp directory to simulate directory input."""
    import shutil
    pdb_dir = tmp_path / "pdbs"
    pdb_dir.mkdir()
    for name in ["discontinuous_structure.pdb", "rebased.pdb"]:
        src = os.path.join("tests/test_data/structure_related", name)
        shutil.copy(src, pdb_dir / name)
    return str(pdb_dir)


def test_directory_mode_print(tmp_path, capsys):
    pdb_dir = _make_test_dir(tmp_path)

    class args:
        structure_file = pdb_dir
        out_file = ""
        header = ""
        threads = 1

    struc_to_seq_main(args)
    captured = capsys.readouterr()
    assert captured.out.count(">") == 2
    assert ">discontinuous_structure\n" in captured.out
    assert ">rebased\n" in captured.out


def test_directory_mode_write(tmp_path):
    pdb_dir = _make_test_dir(tmp_path)
    out = str(tmp_path / "combined.fasta")

    class args:
        structure_file = pdb_dir
        out_file = out
        header = ""
        threads = 1

    struc_to_seq_main(args)

    with open(out) as f:
        content = f.read()
    assert content.count(">") == 2
    assert ">discontinuous_structure\n" in content
    assert ">rebased\n" in content
    seq_lines = [l for l in content.strip().split("\n") if not l.startswith(">")]
    for seq in seq_lines:
        assert seq == EXPECTED_SEQ


def test_directory_mode_overwrites(tmp_path):
    pdb_dir = _make_test_dir(tmp_path)
    out = str(tmp_path / "combined.fasta")

    class args:
        structure_file = pdb_dir
        out_file = out
        header = ""
        threads = 1

    struc_to_seq_main(args)
    struc_to_seq_main(args)

    with open(out) as f:
        content = f.read()
    assert content.count(">") == 2


def test_directory_mode_threads(tmp_path):
    pdb_dir = _make_test_dir(tmp_path)
    out = str(tmp_path / "combined.fasta")

    class args:
        structure_file = pdb_dir
        out_file = out
        header = ""
        threads = 2

    struc_to_seq_main(args)

    with open(out) as f:
        content = f.read()
    assert content.count(">") == 2


def test_directory_mode_failed_file(tmp_path):
    pdb_dir = _make_test_dir(tmp_path)
    bad = os.path.join(pdb_dir, "broken.pdb")
    with open(bad, "w") as f:
        f.write("not a pdb\n")
    out = str(tmp_path / "combined.fasta")

    class args:
        structure_file = pdb_dir
        out_file = out
        header = ""
        threads = 1

    with pytest.raises(ValueError, match="broken.pdb"):
        struc_to_seq_main(args)
