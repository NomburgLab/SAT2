from pathlib import Path
import filecmp
from sat.scripts.struc_get_domains import (
    format_chainsaw_domains,
    parse_domain,
    struc_extract_residues,
    struc_get_domains_main
)

class Test_Parse_Domain():
    def test_parse_single_domain(self):
        domain_boundaries = '1-10'
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        assert parse_domain(domain_boundaries) == expected

class Test_Format_Chainsaw_Domains():

    def test_single_domain_formatting(self):
        domain_dict = {'structure1': {'chopping':"3-8"}}
        expected = list(zip(["3-8"],[[3, 4, 5, 6, 7, 8]]))
        assert list(format_chainsaw_domains('structure1', domain_dict)) == expected

    def test_multiple_domains_formatting(self):
        domain_dict = {'structure1': {'chopping':"3-8,11-15"}}
        expected = list(zip(["3-8", "11-15"] ,[[3, 4, 5, 6, 7, 8], [11, 12, 13, 14, 15]]))
        assert list(format_chainsaw_domains('structure1', domain_dict)) == expected

    def test_single_domain_with_subdomain_formatting(self):
        domain_dict = {'structure1': {'chopping':"3-8_11-15"}}
        expected = list(zip(["3-8_11-15"] ,[[3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15]]))
        assert list(format_chainsaw_domains('structure1', domain_dict)) == expected

    def test_complex_domains_with_subdomains_formatting(self):
        domain_dict = {'structure1': {'chopping':"3-8_11-15,20-25,30-32_40-45"}}
        expected = list(zip( ["3-8_11-15", "20-25", "30-32_40-45"],
                            [[3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15], [20, 21, 22, 23,24,25], [30,31,32,40,41,42,43,44,45]]))
        assert list(format_chainsaw_domains('structure1', domain_dict)) == expected

class Test_Struc_Extract_Residues():

    def test_min_length(self,tmp_path):

        #define input parameters
        pdb_file_path = (
            "tests/test_data/structure_related/get_domains_chainsaw/input_min_length/AB537968__BAJ06111.1__X__00001.pdb"
        )
        zipped_domains_and_residues = zip(
            ["4","4-7","4-8","50-60","7-8_15-16","60-65_70-75_80-85"],
            [[4],
            [4,5,6,7],
            [4,5,6,7,8],
            [50,51,52,53,54,55,56,57,58,59,60],
            [7,8,15,16],
            [60,61,62,63,64,65,70,71,72,73,74,75,80,81,82,83,84,85]]
        )
        min_domain_length = 5
        outfile_dir = tmp_path/"test_output"
        outfile_dir.mkdir()
 
        expected_outfile_dir = Path("tests/test_data/structure_related/get_domains_chainsaw/output_min_length")

        # Running script
        struc_extract_residues(pdb_file_path, zipped_domains_and_residues,  min_domain_length, outfile_dir)

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