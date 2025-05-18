#!/usr/bin/env python3
"""
Convert DALI alignments to pairwise FASTA files.

This script takes a DALI alignment output file and converts each pairwise alignment 
into a separate FASTA file for easier downstream analysis.
"""

import os
import re
import argparse
from typing import Dict, List, Tuple

from .utils.misc import talk_to_me, make_output_dir
from .utils.dali import parse_structure_key, read_alignment_block, segment_alignments, DALI_alignment


def aln_dali_to_pairwise_fastas_main(args):
    """
    Main function to convert DALI alignments to pairwise FASTA files.
    
    Args:
        args: Command-line arguments
    """
    talk_to_me(f"Reading DALI alignment file: {args.aln_file}")
    
    if args.key != "":
        key = parse_structure_key(args.key)
        talk_to_me(f"Using structure key: {args.key}")
    else:
        key = ""
        talk_to_me("No structure key provided")
    
    # Read and segment the alignments
    alignments = read_alignment_block(args.aln_file)
    alignments = segment_alignments(alignments)
    talk_to_me(f"Found {len(alignments)} pairwise alignments")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    talk_to_me(f"Writing pairwise FASTA files to: {args.out_dir}")
    
    # Track stats
    total_alignments = len(alignments)
    skipped_self_alignments = 0
    
    # Process each alignment
    for i, alignment in enumerate(alignments):
        # Parse the alignment
        dali_aln = DALI_alignment(alignment, key)
        dali_aln.parse_alignment()
        
        # Get query and target IDs
        query_id = dali_aln.query
        target_id = dali_aln.target
        
        # Skip self-alignments if requested
        if args.skip_self and query_id == target_id:
            skipped_self_alignments += 1
            continue
        
        # Create output filename
        output_file = os.path.join(args.out_dir, f"{query_id}xxx{target_id}.fasta")
        
        # Write the pairwise FASTA
        with open(output_file, 'w') as f:
            # Write query sequence
            f.write(f">{query_id}\n")
            f.write(f"{dali_aln.aln_qseq}\n")
            
            # Write target sequence
            f.write(f">{target_id}\n")
            f.write(f"{dali_aln.aln_tseq}\n")
        
        # Log progress periodically
        if i < 5 or i == len(alignments) - 1:
            talk_to_me(f"Wrote: {os.path.basename(output_file)}")
        elif i == 5:
            talk_to_me(f"... and {len(alignments) - 6 - skipped_self_alignments} more files")
    
    talk_to_me(f"Completed processing {total_alignments} alignments")
    if args.skip_self and skipped_self_alignments > 0:
        talk_to_me(f"Skipped {skipped_self_alignments} self-alignments")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg) 