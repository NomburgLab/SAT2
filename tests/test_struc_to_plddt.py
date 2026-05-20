import pytest
from sat.scripts.struc_to_plddt import struc_to_plddt_main

SINGLE_PDB = "tests/test_data/structure_related/discontinuous_structure.pdb"
MULTI_PDB_DIR = "tests/test_data/structure_related/get_domains_chainsaw/input_extract_all"


def _make_args(**kwargs):
    class Args:
        pass

    args = Args()
    args.structure_file = ""
    args.out_file = ""
    args.threads = 1
    for k, v in kwargs.items():
        setattr(args, k, v)
    return args


class TestStrucToPlddt:
    def test_single_file(self, capsys):
        args = _make_args(structure_file=SINGLE_PDB)
        struc_to_plddt_main(args)
        output = capsys.readouterr().out.strip().split("\n")[0]
        plddt = float(output)
        assert 0 < plddt <= 100

    def test_directory(self, tmp_path):
        out_file = str(tmp_path / "results.tsv")
        args = _make_args(structure_file=MULTI_PDB_DIR, out_file=out_file)
        struc_to_plddt_main(args)

        with open(out_file) as f:
            lines = f.read().strip().split("\n")
        assert len(lines) == 3
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 2
            assert parts[0].endswith(".pdb")
            plddt = float(parts[1])
            assert 0 < plddt <= 100

    def test_parallel_matches_sequential(self, tmp_path):
        out_seq = str(tmp_path / "seq.tsv")
        args_seq = _make_args(structure_file=MULTI_PDB_DIR, out_file=out_seq, threads=1)
        struc_to_plddt_main(args_seq)

        out_par = str(tmp_path / "par.tsv")
        args_par = _make_args(structure_file=MULTI_PDB_DIR, out_file=out_par, threads=2)
        struc_to_plddt_main(args_par)

        with open(out_seq) as f:
            seq_lines = sorted(f.readlines())
        with open(out_par) as f:
            par_lines = sorted(f.readlines())
        assert seq_lines == par_lines

    def test_output_to_file(self, tmp_path):
        out_file = str(tmp_path / "out.tsv")
        args = _make_args(structure_file=SINGLE_PDB, out_file=out_file)
        struc_to_plddt_main(args)

        with open(out_file) as f:
            content = f.read().strip()
        parts = content.split("\t")
        assert len(parts) == 2
        assert parts[0] == "discontinuous_structure.pdb"
        float(parts[1])
