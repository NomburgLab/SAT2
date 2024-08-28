import pytest
from sat.scripts.aln_trim_target_fasta import (
    select_alignment,
    trim_sequence
)

def test_select_alignment_one_shortest():
    aln_list = [[1,4,4], [2,6,5], [1,2,2]] 
    short_aln_type = True
    expected = [1,2,2]
    assert select_alignment(aln_list, short_aln_type) == expected

def test_select_alignment_one_longest():
    aln_list = [[1,4,4], [2,6,5], [1,2,2], [3,30,28]] 
    short_aln_type = False
    expected = [3,30,28]
    assert select_alignment(aln_list, short_aln_type) == expected

def test_select_alignment_multiple_shortest():
    aln_list = [[2,4,3], [3,5,3], [2,5,4], [3,5,3]]
    short_aln_type = True
    expected = [2,4,3]
    assert select_alignment(aln_list, short_aln_type) == expected

def test_select_alignment_mutiple_longest():
    aln_list = [[3,30,28], [2,6,5], [2,29,28], [1,2,2], [3,30,28]] 
    short_aln_type = False
    expected = [3,30,28]
    assert select_alignment(aln_list, short_aln_type) == expected

def test_trim_sequence_empty_aln_list():
    chosen_aln_list = [1,2]
    sequence = '123456789'
    with pytest.raises(ValueError):
        trim_sequence(chosen_aln_list, sequence )

def test_trim_sequence_correct_aln_length():
    sequence = '123456789'
    chosen_aln_list = [3, 7, 5]
    expected = '34567'
    assert expected == trim_sequence(chosen_aln_list, sequence)

def test_trim_sequence_wrong_aln_length():
    sequence = '123456789'
    chosen_aln_list = [3, 7, 4]
    with pytest.raises(ValueError):
        trim_sequence(chosen_aln_list, sequence )