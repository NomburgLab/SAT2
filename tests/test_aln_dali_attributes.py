import pytest
from sat.scripts.aln_dali_attributes import (
    DALI_alignment_attributes
)

def test_find_aln_position():
    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  ---ABC-D--EF-G-Fy   9",  
    "ident     ||| |     |       ",   
    "Sbjct  --dABC-D--XY-G--y   9",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")
    aln_obj.find_aln_position()

    assert aln_obj.aln_start == 3
    assert aln_obj.aln_end == 13

def test_find_aln_position_no_aln():
    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  ---a-----------fy   3",  
    "ident                       ",   
    "Sbjct  --d-bd-d--yx-g--y   8",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")

    with pytest.raises(ValueError):
        aln_obj.find_aln_position()

def test_find_seq_position():

    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  ---ABC-D--EF-G-fy   9",  
    "ident     ||| |     |       ",   
    "Sbjct  --dABC-D--XY-G---   8",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")
    seq_position = aln_obj.find_seq_position(7, aln_obj.aln_qseq)

    assert seq_position == 4

def test_find_seq_position_out_of_range():
    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  ---ABC-D--EF-G-fy   9",  
    "ident     ||| |     |       ",   
    "Sbjct  --dABC-D--XY-G---   8",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")
    with pytest.raises(ValueError):
        aln_obj.find_seq_position(18, aln_obj.aln_qseq)

def test_calculate_pident():
    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  -------D--EF-G-FY   6",  
    "ident         |     |       ",   
    "Sbjct  --dabc-D--XY-G-cd   10",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")
    aln_obj.pident = aln_obj.calculate_pident()
    assert aln_obj.pident == 50

def test_calculate_pident_no_aligned_residues():
    aln= [
    "No 1: Query=ABC Sbjct=DABC Z-score=9.3",
    "DSSP  -----L--L---L-L--",   
    "Query  ------c--d--v--fy   5",  
    "ident                       ",   
    "Sbjct  --dabc-d--yx-g-fy   10",  
    "DSSP  ----L--L---L-L---"     
    ]

    aln_obj = DALI_alignment_attributes(aln, key ="")
    with pytest.raises(ValueError):
        aln_obj.calculate_pident()