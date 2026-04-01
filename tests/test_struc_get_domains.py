from pathlib import Path
import filecmp
import pytest
from sat.scripts.struc_get_domains import (
    parse_domain,
    struc_get_domains_main
)

# Shared paths
INPUT_EXTRACT_ALL = Path("tests/test_data/structure_related/get_domains_chainsaw/input_extract_all")
OUTPUT_EXTRACT_ALL = Path("tests/test_data/structure_related/get_domains_chainsaw/output_extract_all")
INPUT_MERIZO = Path("tests/test_data/structure_related/get_domains_chainsaw/input_merizo")
OUTPUT_MERIZO = Path("tests/test_data/structure_related/get_domains_chainsaw/output_merizo")

PDB_MN539721 = str(INPUT_EXTRACT_ALL / "MN539721__QGH71255.1__X__00025.pdb")
PDB_MN876845 = str(INPUT_EXTRACT_ALL / "MN876845__QJF12414.1__PSV2-gp02__00002.pdb")
PDB_MT047590 = str(INPUT_EXTRACT_ALL / "MT047590__QIM61606.1__X__00001.pdb")


def _make_args(**kwargs):
    """Create a simple args namespace from keyword arguments."""
    class args:
        pass
    for k, v in kwargs.items():
        setattr(args, k, v)
    return args


def _compare_outputs(observed_dir, expected_dir):
    """Assert that every file in observed_dir matches the same-named file in expected_dir."""
    observed_files = list(Path(observed_dir).iterdir())
    assert observed_files, "No output files were produced"
    for observed_outfile in observed_files:
        expected_file = Path(expected_dir) / observed_outfile.name
        assert expected_file.exists(), f"Expected file not found: {expected_file}"
        assert filecmp.cmp(expected_file, observed_outfile, shallow=False), \
            f"Output mismatch: {observed_outfile.name}"


class Test_Parse_Domain():
    def test_parse_single_domain(self):
        domain_boundary = '1-10'
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        assert parse_domain(domain_boundary) == expected

    def test_single_domain_with_subdomain_formatting(self):
        domain_boundary = "3-8_11-15"
        expected = [3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15]
        assert parse_domain(domain_boundary) ==  expected

    def test_multiple_domains(self):
        domain_boundary= "3-8,11-15"
        with pytest.raises(ValueError):
            parse_domain(domain_boundary)

    def test_single_residue_domain(self):
        domain_boundary= "3"
        with pytest.raises(ValueError):
            parse_domain(domain_boundary)

    def test_empty_domain(self):
        domain_boundary= " "
        with pytest.raises(ValueError):
            parse_domain(domain_boundary)

    def test_no_hyphen_domain(self):
        domain_boundary = "10_20"
        with pytest.raises(ValueError):
            parse_domain(domain_boundary)


class Test_Struc_Get_Domains_Main():

    def test_min_length(self, tmp_path):
        args = _make_args(
            structure_file_path="tests/test_data/structure_related/get_domains_chainsaw/input_min_length/AB537968__BAJ06111.1__X__00001.pdb",
            domain_file_path="tests/test_data/structure_related/get_domains_chainsaw/input_min_length/chainsaw_file_min_length_test.txt",
            colnames="",
            domain_column="chopping",
            id_column=None,
            min_domain_length=5,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        struc_get_domains_main(args)

        _compare_outputs(
            args.outfile_dir,
            "tests/test_data/structure_related/get_domains_chainsaw/output_min_length",
        )

    def test_extract_all_tsv(self, tmp_path):
        """Standard Chainsaw TSV with header. NULL chopping → __DUNK."""
        args = _make_args(
            domain_file_path=str(INPUT_EXTRACT_ALL / "chainsaw_file_extract_all_test.txt"),
            colnames="",
            domain_column="chopping",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        for pdb in INPUT_EXTRACT_ALL.iterdir():
            if pdb.suffix == '.pdb':
                args.structure_file_path = str(pdb)
                struc_get_domains_main(args)

        _compare_outputs(args.outfile_dir, OUTPUT_EXTRACT_ALL)

    def test_no_header_tsv(self, tmp_path):
        """TSV without a header row: same data as extract_all, tested on MN876845."""
        args = _make_args(
            structure_file_path=PDB_MN876845,
            domain_file_path=str(INPUT_EXTRACT_ALL / "chainsaw_file_extract_all_no_header.txt"),
            colnames="chain_id,sequence_md5,nres,ndom,chopping,confidence,time_sec",
            domain_column="chopping",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        struc_get_domains_main(args)
        _compare_outputs(args.outfile_dir, OUTPUT_EXTRACT_ALL)

    def test_missing_structure_raises(self, tmp_path):
        """Structure not in domain file → ValueError."""
        args = _make_args(
            structure_file_path=PDB_MN539721,
            domain_file_path=str(INPUT_EXTRACT_ALL / "chainsaw_file_missing_structure.txt"),
            colnames="",
            domain_column="chopping",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        with pytest.raises(ValueError, match="not found in the domain file"):
            struc_get_domains_main(args)

    def test_merizo_format(self, tmp_path):
        """Merizo headerless TSV: IDs include .pdb extension, boundaries in last column."""
        args = _make_args(
            domain_file_path=str(INPUT_MERIZO / "merizo_test.txt"),
            colnames="id,nres,nres_resolved,nres_unresolved,ndom,confidence,mean_plddt,domains",
            domain_column="domains",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        for pdb in INPUT_MERIZO.glob("*.pdb"):
            args.structure_file_path = str(pdb)
            struc_get_domains_main(args)

        _compare_outputs(args.outfile_dir, OUTPUT_MERIZO)

    def test_full_chain_domain_produces_dfull(self, tmp_path):
        """Explicit boundary spanning all PDB residues → output __DFULL."""
        args = _make_args(
            structure_file_path=PDB_MN539721,
            domain_file_path=str(INPUT_EXTRACT_ALL / "chainsaw_file_full_chain.txt"),
            colnames="",
            domain_column="chopping",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        struc_get_domains_main(args)

        output_files = list(args.outfile_dir.iterdir())
        assert len(output_files) == 1
        assert output_files[0].name == "MN539721__QGH71255.1__X__00025__DFULL.pdb"
        expected = OUTPUT_EXTRACT_ALL / "MN539721__QGH71255.1__X__00025__DFULL.pdb"
        assert filecmp.cmp(expected, output_files[0], shallow=False)

    def test_directory_mode(self, tmp_path):
        """Passing a directory of PDB files processes all of them."""
        args = _make_args(
            structure_file_path=str(INPUT_EXTRACT_ALL),
            domain_file_path=str(INPUT_EXTRACT_ALL / "chainsaw_file_extract_all_test.txt"),
            colnames="",
            domain_column="chopping",
            id_column=None,
            min_domain_length=1,
            outfile_dir=tmp_path / "test_output",
        )
        args.outfile_dir.mkdir()

        struc_get_domains_main(args)

        _compare_outputs(args.outfile_dir, OUTPUT_EXTRACT_ALL)
