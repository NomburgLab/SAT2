from pathlib import Path
import filecmp
import pytest
from sat.scripts.struc_get_domains import (
    parse_domain,
    struc_extract_residues,
    struc_get_domains_main
)

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


class Test_Struc_Extract_Residues():

    def test_min_length(self,tmp_path):

        #define input parameters
        pdb_file_path = (
            "tests/test_data/structure_related/get_domains_chainsaw/input_min_length/AB537968__BAJ06111.1__X__00001.pdb"
        )
        domain_dict = {"AB537968__BAJ06111.1__X__00001":
                       {'chopping':'4-7,4-8,50-60,7-8_15-16,60-65_70-75_80-85', 'nres':94}}
        min_domain_length = 5
        outfile_dir = tmp_path/"test_output"
        outfile_dir.mkdir()
 
        expected_outfile_dir = Path("tests/test_data/structure_related/get_domains_chainsaw/output_min_length")

        # Running script
        struc_extract_residues(pdb_file_path, domain_dict,  min_domain_length, outfile_dir)

        # Compare observed and expected
        for observed_outfile in  outfile_dir.iterdir():
            expected_file = expected_outfile_dir/observed_outfile.name
            assert filecmp.cmp(expected_file,observed_outfile, shallow=False)
    
    def test_struc_get_domains_main(self, tmp_path):

        #define input parameters
        class args:
            pass

        args.chainsaw_file_path = (
            "tests/test_data/structure_related/get_domains_chainsaw/input_main/chainsaw_file_test.txt"
        )
        args.min_domain_length = 1
        args.outfile_dir = tmp_path/"test_output"
        args.outfile_dir.mkdir()
        expected_outfile_dir = Path("tests/test_data/structure_related/get_domains_chainsaw/output_main")

        # Running script
        pdb_file_folder = Path("tests/test_data/structure_related/get_domains_chainsaw/input_main/")
        for file in pdb_file_folder.iterdir():
            if file.suffix == '.pdb':
                args.structure_file_path = str(file)
                struc_get_domains_main(args)

        # Compare observed and expected
        for observed_outfile in args.outfile_dir.iterdir():
            expected_file = expected_outfile_dir/observed_outfile.name
            assert filecmp.cmp(expected_file,observed_outfile, shallow=False)