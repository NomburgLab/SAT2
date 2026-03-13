import os

import pytest

from sat.scripts.struc_clean import struc_clean_main # type: ignore


# ------------------------------------------------------------------------------------ #
# Whole-script tests
# ------------------------------------------------------------------------------------ #
def test_struc_clean_main(tmp_path):

    # Define inputs
    class args:
        pass

    args.infile  = "tests/test_data/structure_related/clean/1NWR.pdb"
    args.outfile = f"{tmp_path}/1NWR_clean.pdb"

    # Run the script
    struc_clean_main(args)

    # Validate that the output file was created
    assert os.path.exists(args.outfile), "Output PDB file was not created."

    with open(args.outfile) as f:
        actual_lines = f.readlines()

    # Validate no non-protein residues slipped through
    from sat.scripts.struc_clean import STANDARD_AA
    for line in actual_lines:
        record = line[:6].strip()
        if record in ("ATOM", "HETATM", "ANISOU"):
            resname = line[17:20].strip()
            assert resname in STANDARD_AA, (
                f"Non-protein residue found in output: {resname!r} — line: {line.rstrip()}"
            )

    # Validate the expected set of residues is present (same chain/resseq as reference)
    def extract_residues(path):
        residues = set()
        with open(path) as f:
            for line in f:
                if line[:6].strip() in ("ATOM", "HETATM"):
                    chain  = line[21]
                    resseq = line[22:26].strip()
                    resname = line[17:20].strip()
                    residues.add((chain, resseq, resname))
        return residues

    expected_residues = extract_residues("tests/test_data/structure_related/clean/1NWR_noligands.pdb")
    actual_residues   = extract_residues(args.outfile)
    assert actual_residues == expected_residues, (
        f"Residue mismatch.\n"
        f"  Only in actual:   {actual_residues - expected_residues}\n"
        f"  Only in expected: {expected_residues - actual_residues}"
    )


# ------------------------------------------------------------------------------------ #
# Edge-case / error-handling tests
# ------------------------------------------------------------------------------------ #
def test_struc_clean_bad_extension(tmp_path):
    """Should raise ValueError when infile does not end with .pdb."""

    class args:
        pass

    args.infile  = "tests/test_data/structure_related/clean/1NWR.cif"
    args.outfile = f"{tmp_path}/out.pdb"

    with pytest.raises(ValueError, match=r"\.pdb"):
        struc_clean_main(args)


def test_struc_clean_missing_file(tmp_path):
    """Should raise ValueError when infile does not exist."""

    class args:
        pass

    args.infile  = "tests/test_data/structure_related/clean/DOES_NOT_EXIST.pdb"
    args.outfile = f"{tmp_path}/out.pdb"

    with pytest.raises(ValueError, match="Cannot detect PDB file"):
        struc_clean_main(args)