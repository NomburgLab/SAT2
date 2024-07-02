from sat.scripts.struc_get_domains import (
    format_chainsaw_domains,
    struc_extract_residues
)

from unittest.mock import patch, MagicMock, call


def test_single_domain_formatting():
    domain_dict = {'structure1': "3-8"}
    expected = [[3, 4, 5, 6, 7, 8]]
    assert format_chainsaw_domains('structure1', domain_dict) == expected

def test_multiple_domains_formatting():
    domain_dict = {'structure1': "3-8, 11-15"}
    expected = [[3, 4, 5, 6, 7, 8], [11, 12, 13, 14, 15]]
    assert format_chainsaw_domains('structure1', domain_dict) == expected

def test_single_domain_with_subdomain_formatting():
    domain_dict = {'structure1': "3-8_11-15"}
    expected = [[3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15]]
    assert format_chainsaw_domains('structure1', domain_dict) == expected

def test_complex_domains_with_subdomains_formatting():
    domain_dict = {'structure1': "3-8_11-15, 20-25, 30-32_40-45"}
    expected = [[3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15], [20, 21, 22, 23,24,25], [30,31,32,40,41,42,43,44,45]]
    assert format_chainsaw_domains('structure1', domain_dict) == expected

@patch('struc_get_domains.pdb_to_structure_object')
@patch('struc_get_domains.get_filename')
@patch('struc_get_domains.write_structure_subset')
def test_domain_length_less_than_min(mock_write_structure_subset, mock_get_filename, mock_pdb_to_structure_object):

    mock_structure = MagicMock()
    mock_pdb_to_structure_object.return_value = mock_structure
    mock_get_filename.return_value = 'test_structure'
    
    pdb_file_path = 'test.pdb'
    domain_residues_list = [
        [1],
        [1,2],
        [1,2,3],
        [1,2,3,4],
        [2,3,7,8,9]
    ]
    min_domain_length = 3
    output_dir = 'output'

    struc_extract_residues(pdb_file_path, domain_residues_list, min_domain_length, output_dir)

    expected_calls = [
        call(mock_structure, residues_to_keep=[1,2,3], outfile='output/test_structure_domain_1_3.pdb'),
        call(mock_structure, residues_to_keep=[1,2,3,4], outfile='output/test_structure_domain_1_4.pdb'),
        call(mock_structure, residues_to_keep=[2,3,7,8,9], outfile='output/test_structure_domain_2_9.pdb')
    ]
    mock_write_structure_subset.assert_has_calls(expected_calls, any_order=False)
    assert mock_write_structure_subset.call_count == 3

    