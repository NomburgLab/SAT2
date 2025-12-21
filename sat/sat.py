#!/usr/bin/env python

import argparse


def arg_str2bool(v):
    """
    For use as an argparse argument type. Makes it easy to use boolean flags.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def main():

    # Top-level parser
    parser = argparse.ArgumentParser(
        description=(
            "SAT - Structural Analysis Toolkit. A python package for manipulating "
            "predicted structures and structural alignments."
        ),
        usage="""sat.py""",
    )
    subparsers = parser.add_subparsers(
        title="Subcommands",
        required=True,
    )

    # -------------------------------------------------------------------------------- #
    # Parser for DALI_alignment_attributes subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_dali_alignment_attributes = subparsers.add_parser(
        "aln_dali_alignment_attributes",
        help=(
            """
            This subcommands takes in a DALI alignment file (which must have an alignments field) 
            and a key to generate a csv file of aligned targets and queries and their attributes, 
            such as qstart, qend, tstart, tend, etc.
            """
        ),
    )
    parser_aln_dali_alignment_attributes.add_argument(
        "-i",
        "--aln_file",
        type=str,
        required=True,
        help="""
        Path to a DALI file that contains the 'alignments' output field. There can be
        other output fields as well, but only the alignments will be parsed here.
        """,
    )
    parser_aln_dali_alignment_attributes.add_argument(
        "-k",
        "--key",
        type=str,
        required=False,
        default="",
        help="""
        Path to a tab-delimited file of format structure_name,,ID, where
        the ID is a 4-digit identifier used during the DALI alignment. This lets you
        convert the identifiers back.
        """,
    )
    parser_aln_dali_alignment_attributes.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="""
        Path to the output file. This file will have two columns: target and target_id.
        The target_id is the DALI ID, while the target will be present if a key was
        used.
        """,
    )
    parser_aln_dali_alignment_attributes.set_defaults(
        func=call_aln_dali_alignment_attributes_main
    )

    # -------------------------------------------------------------------------------- #
    # Parser for aln_add_clusters subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_add_clusters = subparsers.add_parser(
        "aln_add_clusters",
        help=(
            """
            This subcommand incorporates clustering information from Foldseek cluster
            into the foldseek alignment tabular output file. Notably, this script
            generates "super clusters" based on the foldseek cluster file. Here,
            we determine what other clusters each cluster is linked to - two clusters
            are linked if at least linkage_threshold fraction of at least one of the
            cluster's members have alignments to members of the other cluster. Then,
            "super clusters" are generated from clusters whom are all linked to one
            another. Thus, the cluster outputs from this script reflect the super
            clusters. The cluster representative from the foldseek cluster with the most
            members is then used as the super cluster representative.

            This script makes two output files:
            1) An output alignment file with all non-redundant alignments.
            2) An output alignment file with only the alignments with the cluster
               representative as the query.

            Notably, the new columns "cluster_ID", "cluster_count", and
            "cluster_rep" are added to the output files.
            """
        ),
    )
    parser_aln_add_clusters.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek alignment file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-c",
        "--cluster_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the foldseek cluster tsv file. Typically this will be generated using
        cluster mode 0.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-r",
        "--rep_out",
        type=str,
        required=True,
        help="""
        Path to the output file. Here, for each cluster only alignments with the
        foldseek cluster representative as the query will be output. This is nice
        because the output file size is much smaller than the non-redundant file.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-n",
        "--nr_out",
        type=str,
        required=True,
        help="""
        Path to the output file. This file contains all alignments, although if there
        were self alignments in the input those are removed.
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        A comma-delimited string of the fields in the input foldseek alignment file.
        Make sure to wrap in quotes!

        Default: ''
        """,
    )
    parser_aln_add_clusters.add_argument(
        "-l",
        "--linkage_threshold",
        type=float,
        required=False,
        default=0.5,
        help="""
        [Default: 0.5] The fraction of members of at least one cluster whom must have
        alignments to members of another cluster for those clusters to be considered
        linked.
        """,
    )
    parser_aln_add_clusters.set_defaults(func=call_aln_add_clusters)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_add_uniprot subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_add_uniprot = subparsers.add_parser(
        "aln_add_uniprot",
        help=(
            """
            This subcommand takes in uniprot information generated by aln_query_uniprot
            and adds it to a foldseek alignment.
            """
        ),
    )
    parser_aln_add_uniprot.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. It is OK if the first row is the header, as long as
        the first column is 'query'.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-u",
        "--uniprot_information",
        type=str,
        required=True,
        help="""
        Path to the uniprot information file generated by aln_query_uniprot.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="""
        Path to the output file.
        """,
    )
    parser_aln_add_uniprot.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the input file. Note that if there are labeled
        column headers in the input alignment file, they can be automatically parsed
        and you can leave the alignment_fields command empty.
        """,
    )
    parser_aln_add_uniprot.set_defaults(func=call_aln_add_uniprot_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_connection_map subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_connection_map = subparsers.add_parser(
        "aln_connection_map",
        help=(
            """
            This subcommand takes in a cluster file that has taxonomy information
            (critically - family) and determines, for each pair of families,
            how many clusters exist in which both families have a member.
            """
        ),
    )
    parser_aln_connection_map.add_argument(
        "-c",
        "--cluster_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to cluster file. Taxonomy information is required! There should be a
        species and a family column, at least.s
        """,
    )
    parser_aln_connection_map.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        default="",
        help="""
        Path to output connection map. There are the following columns:
         ["f1", "f2", "count", "f1_total", "f2_total", "jaccard"]
        1. Family 1
        2. Family 2
        3. Count - Number of clusters in which a member of family 1 and family 2 are
            both present.
        4. f1_total: Total number of clusters with a member from family 1
        5. f2_total: Total number of clusters with a member from family 2
        6. Jaccard index: Count/(f1_total + f2_total - Count)
        """,
    )
    parser_aln_connection_map.add_argument(
        "-F",
        "--cluster_file_fields",
        type=str,
        required=False,
        default="",
        help="""
        [Default: '']
        Comma-delimited string indicating the columns present in the cluster_file.
        """,
    )
    parser_aln_connection_map.set_defaults(func=call_aln_connection_map)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_dali_motif_finder subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_dali_motif_finder = subparsers.add_parser(
        "aln_dali_motif_finder",
        help=(
            """
            This subcommand parses a DALI alignments field (specified by the
            'alignments' output format) and determines, for each target, if a specified
            motif or series of motifs is present at specified locations.

            The output of this script is a a file with two columns - target and 
            target_id - for those targets that contain all indicated motifs. Note that
            adding a DALI key is optional, in which case 'target' will be blank.
             
            See below for detailed instructions about how to input motifs.
            """
        ),
    )
    parser_aln_dali_motif_finder.add_argument(
        "-i",
        "--aln_file",
        type=str,
        required=True,
        help="""
        Path to a DALI file that contains the 'alignments' output field. There can be
        other output fields as well, but only the alignments will be parsed here.
        """,
    )
    parser_aln_dali_motif_finder.add_argument(
        "-m",
        "--motif_list",
        type=str,
        required=True,
        help="""
        See the documentation on github for a detailed description.
        """,
    )
    parser_aln_dali_motif_finder.add_argument(
        "-k",
        "--key",
        type=str,
        required=False,
        default="",
        help="""
        Path to a tab-delimited file of format structure_name,,ID, where
        the ID is a 4-digit identifier used during the DALI alignment. This lets you
        convert the identifiers back.
        """,
    )
    parser_aln_dali_motif_finder.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="""
        Path to the output file. This file will have two columns: target and target_id.
        The target_id is the DALI ID, while the target will be present if a key was
        used.
        """,
    )
    parser_aln_dali_motif_finder.set_defaults(func=call_aln_dali_motif_finder_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_dali_to_pairwise_fastas subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_dali_to_pairwise_fastas = subparsers.add_parser(
        "aln_dali_to_pairwise_fastas",
        help="""
        This subcommand converts the alignment fields of a DALI output file to pairwise
        fasta files, one per alignment. The usecase here is typically when you are
        trying to generate a structure-guided MSA. This means that, during the DALI 
        run, the 'alignments' output field must have been included.
        
        For each pairwise alignment in the DALI output, a separate FASTA file is created
        with the naming pattern [query_id]xxx[target_id].fasta.aln in the specified output directory.
        Each FASTA file contains two sequences:
        1. The query sequence from the alignment
        2. The target sequence from the alignment
        
        If a structure key file is provided, DALI IDs will be converted to the corresponding
        structure names found in the key file.
        """,
    )
    parser_aln_dali_to_pairwise_fastas.add_argument(
        "-a",
        "--aln_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input DALI alignment file.
        """,
    )
    parser_aln_dali_to_pairwise_fastas.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=True,
        default="",
        help="""
        Directory to write the output FASTA files to.
        """,
    )
    parser_aln_dali_to_pairwise_fastas.add_argument(
        "-k",
        "--key",
        type=str,
        required=False,
        default="",
        help="""
        Path to a structure key file to convert DALI IDs to structure names.
        """,
    )
    parser_aln_dali_to_pairwise_fastas.add_argument(
        "-s",
        "--skip_self",
        action="store_true",
        help="""
        Skip self-alignments (where query and target are the same).
        By default, all alignments including self-alignments are processed.
        """,
    )

    parser_aln_dali_to_pairwise_fastas.set_defaults(func=call_aln_dali_to_pairwise_fastas_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_filter subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_filter = subparsers.add_parser(
        "aln_filter",
        help=(
            """
            This subcommand filters a foldseek alignment file to keep only those
            alignments with a value below/above the specified value in a field
            (alntmscore is a common one). It also only outputs a maximum of N alignments
            for each query.
            """
        ),
    )
    parser_aln_filter.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. It is OK if the first row is the header, as long as
        the first column is 'query'.
        """,
    )
    parser_aln_filter.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output alignment file.
        """,
    )
    parser_aln_filter.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of alignment fields. Leave blank if there is a header in
        the input file.
        """,
    )
    parser_aln_filter.add_argument(
        "-x",
        "--filter_field",
        type=str,
        required=False,
        default="alntmscore",
        help="""
        Default: 'alntmscore'
        String indicating the alignment field which contains the value to be filtered/
        sorted with.
        """,
    )
    parser_aln_filter.add_argument(
        "-m",
        "--min_val_filter_field",
        type=float,
        required=False,
        default=0.4,
        help="""
        Default: 0.4
        An alignment must have at least this minimum value in the filter_field for the
        alignment to be output.
        """,
    )
    parser_aln_filter.add_argument(
        "-M",
        "--max_val_filter_field",
        type=float,
        required=False,
        default=1,
        help="""
        Default: 1
        An alignment must have this value or less in the filter_field for the
        alignment to be output.
        """,
    )
    parser_aln_filter.set_defaults(func=call_aln_filter_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_cluster_greedy subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_cluster_greedy = subparsers.add_parser(
        "aln_cluster_greedy",
        help=(
            """
            This subcommand performs greedy set cover clustering, similar to
            foldseek/mmseqs cluster mode 0.

            The algorithm:
            1. Sorts members by number of alignments (descending)
            2. Picks the top member as cluster rep, assigns all its alignments to its cluster
            3. Removes assigned members from consideration
            4. Repeats until all members are clustered
            5. Performs a reassignment step where members are moved to a different
               cluster if they have a better alignment score to that cluster's rep.

            Unlike connected-component clustering, this method does NOT transitively
            connect members. If A aligns to B and B aligns to C, but A doesn't align
            to C, then A and C may end up in different clusters.

            The output is a foldseek-style cluster file with two columns:
            cluster_rep and cluster_member.
            """
        ),
    )
    parser_aln_cluster_greedy.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. If the first row starts with 'query', it will
        be automatically used as the header. Otherwise, provide column names via
        --colnames.
        """,
    )
    parser_aln_cluster_greedy.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output cluster file. The file will have two columns:
        cluster_rep and cluster_member.
        """,
    )
    parser_aln_cluster_greedy.add_argument(
        "-c",
        "--colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of column names in the alignment file. If not
        provided, the first line of the file must start with 'query' and will
        be used as the header.
        """,
    )
    parser_aln_cluster_greedy.add_argument(
        "-s",
        "--score_column",
        type=str,
        required=False,
        default="alntmscore",
        help="""
        [Default: 'alntmscore']
        Name of the column containing alignment scores. Higher scores indicate
        better alignments. Used for the reassignment step.
        """,
    )
    parser_aln_cluster_greedy.add_argument(
        "-A",
        "--all_inputs",
        type=str,
        required=False,
        default="",
        help="""
        [Default: '']
        Path to a file that lists, one per line, all members that were initially
        input into the alignment algorithm. Members not in the alignment file
        will be added as single-member clusters.
        """,
    )
    parser_aln_cluster_greedy.add_argument(
        "--no_reassign",
        action="store_true",
        help="""
        Skip the reassignment step. By default, after initial clustering, members
        are reassigned to a different cluster rep if they have a better alignment
        score to that rep.
        """,
    )
    parser_aln_cluster_greedy.set_defaults(func=call_aln_cluster_greedy_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_cigar_to_cov subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_cigar_to_cov = subparsers.add_parser(
        "aln_cigar_to_cov",
        help=(
            """
            This subcommand adds CIGAR-derived coverage columns to an alignment file.
            It parses the CIGAR string to count the number of M (match) operations,
            then calculates query coverage (cigar_qcov = #Ms/qlen) and target coverage
            (cigar_tcov = #Ms/tlen). These columns are inserted after the 'cigar'
            column in the output.

            The input alignment file must have 'qlen', 'tlen', and 'cigar' columns.
            """
        ),
    )
    parser_aln_cigar_to_cov.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. If the first row starts with 'query', it will
        be automatically used as the header. Otherwise, provide column names via
        --colnames.
        """,
    )
    parser_aln_cigar_to_cov.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output alignment file.
        """,
    )
    parser_aln_cigar_to_cov.add_argument(
        "-c",
        "--colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of column names in the alignment file. If not
        provided, the first line of the file must start with 'query' and will
        be used as the header. The column names must include 'qlen', 'tlen',
        and 'cigar'.
        """,
    )
    parser_aln_cigar_to_cov.set_defaults(func=call_aln_cigar_to_cov_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_merge subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_merge = subparsers.add_parser(
        "aln_merge",
        help=(
            """
            This subcommand merges to foldseek alignments. Columns that are present
            in one but not present in the other will be carried over.
            """
        ),
    )
    parser_aln_merge.add_argument(
        "-1",
        "--aln1",
        type=str,
        required=True,
        help="""
        Path the first alignment file.
        """,
    )
    parser_aln_merge.add_argument(
        "-2",
        "--aln2",
        type=str,
        required=True,
        help="""
        Path the second alignment file.
        """,
    )
    parser_aln_merge.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""
        Path the output file.
        """,
    )
    parser_aln_merge.add_argument(
        "-f",
        "--aln1_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the first alignment file. Can leave blank if the
        alignment file has headers.
        """,
    )
    parser_aln_merge.add_argument(
        "-F",
        "--aln2_fields",
        type=str,
        required=False,
        default="",
        help="""
        Alignment fields present in the second alignment file. Can leave blank if the
        alignment file has headers.
        """,
    )
    parser_aln_merge.add_argument(
        "-s",
        "--aln1_source",
        type=str,
        required=False,
        default="",
        help="""
        String indicating the source of alignment 1. This value will be added to each
        alignment that originated from alignment 1. Leave this blank if alignment 1
        already has a source column or you don't care about making a source column.
        """,
    )
    parser_aln_merge.add_argument(
        "-S",
        "--aln2_source",
        type=str,
        required=False,
        default="",
        help="""
        String indicating the source of alignment 2. This value will be added to each
        alignment that originated from alignment 2. Leave this blank if alignment 2
        already has a source column or you don't care about making a source column.
        """,
    )

    parser_aln_merge.set_defaults(func=call_aln_merge_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_merge_clusters subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_merge_clusters = subparsers.add_parser(
        "aln_merge_clusters",
        help=(
            """
            This subcommand merges two cluster files by adding a higher-level clustering
            to a nested cluster file. File1 contains the higher-level clustering
            (e.g., lol_rep -> struc_rep). File2 contains the lower-level/nested
            clustering (e.g., struc_rep -> seq_rep -> member). The second column of
            file1 is joined to the first column of file2. File suffixes (e.g., .pdb,
            .fasta) are automatically removed from all values. 

            If a value in file2's first column is not found in file1's second column,
            it will be used as the higher-level cluster representative. This would be the case,
            for example, if you generated the higher level cluster file from an alignment
            and the protein didn't have any alignments.
            """
        ),
    )
    parser_aln_merge_clusters.add_argument(
        "-1",
        "--file1",
        type=str,
        required=True,
        help="""
        Path to the higher-level cluster file. This file's second column will be
        joined to file2's first column. Can be a simple 2-column file (rep, member)
        or a nested file with more columns.
        """,
    )
    parser_aln_merge_clusters.add_argument(
        "-2",
        "--file2",
        type=str,
        required=True,
        help="""
        Path to the lower-level/nested cluster file. This file's first column will
        be joined to file1's second column. Can be a simple 2-column file or a
        nested file with more columns.
        """,
    )
    parser_aln_merge_clusters.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output merged cluster file.
        """,
    )
    parser_aln_merge_clusters.add_argument(
        "--file1_colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited list of column names for file1. If not provided, the first
        line of file1 will be used as the header.
        """,
    )
    parser_aln_merge_clusters.add_argument(
        "--file2_colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited list of column names for file2. If not provided, the first
        line of file2 will be used as the header.
        """,
    )
    parser_aln_merge_clusters.set_defaults(func=call_aln_merge_clusters_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_cluster_connected_component subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_cluster_connected_component = subparsers.add_parser(
        "aln_cluster_connected_component",
        help=(
            """
            This subcommand generates connected component clusters from an alignment
            file. All query-target pairs that are connected (directly or transitively)
            will be placed into the same cluster. This is similar to foldseek/mmseqs
            cluster mode 1 (connected-component clustering).

            The output is a foldseek-style cluster file with two columns:
            cluster_rep and cluster_member. The cluster representative is selected
            randomly from each cluster.
            """
        ),
    )
    parser_aln_cluster_connected_component.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the alignment file. If the first row starts with 'query', it will
        be automatically used as the header. Otherwise, provide column names via
        --colnames.
        """,
    )
    parser_aln_cluster_connected_component.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output cluster file. The file will have two columns:
        cluster_rep and cluster_member.
        """,
    )
    parser_aln_cluster_connected_component.add_argument(
        "-c",
        "--colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of column names in the alignment file. If not
        provided, the first line of the file must start with 'query' and will
        be used as the header. The column names must include 'query' and 'target'.
        """,
    )
    parser_aln_cluster_connected_component.add_argument(
        "-A",
        "--all_inputs",
        type=str,
        required=False,
        default="",
        help="""
        [Default: '']
        Path to a file that lists, one per line, all members that were initially
        input into the alignment algorithm. Members who did not make an alignment
        (e.g., were filtered out) are not in the alignment file and thus wouldn't
        end up in the cluster file output. By specifying those members here, they
        will be output as single-member clusters.

        The input file should indicate the basename of each structure (if using
        foldseek) or sequence (if using mmseqs) that was initially input into the
        alignment, with one member per line. If the alignment was from foldseek,
        the alignment file queries/targets likely end in .pdb - so the entries in
        this file must also end in .pdb.
        """,
    )
    parser_aln_cluster_connected_component.add_argument(
        "-d",
        "--diagnose",
        action="store_true",
        help="""
        Run diagnostic analysis on the alignment file. This will print information
        about what's connecting clusters together, including the most connected
        members (potential "hub" proteins) and examples of non-self alignments.
        Useful for debugging why you might get one giant cluster.
        """,
    )
    parser_aln_cluster_connected_component.set_defaults(func=call_aln_cluster_connected_component_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_parse_dali subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_parse_dali = subparsers.add_parser(
        "aln_parse_dali",
        help=(
            """
            This subcommand reads in a DALI alignment output file and formats it as
            a tab-delimited file. This script will written to the specified output file.
            There is also functionality to filter the alignments by zscore, alnlen,
            coverage, or rmsd.

            There are two main inputs:
            1) alignment_file: This is the DALI alignment file. Notably, the first
                out put field MUST BE the 'summary' and the second output MUST BE
                'equivalences'.
            2) structure_key: DALI only processes files that have a 4-digit identifier.
                The structure key must be of format structure[delimiter]identifier, and
                lets you convert the identifiers back to the actual structure name.
                Note that the structure_key identifiers should not have the DALI
                segment (e.g. A, B, C...) at the end - this will be taken care of.

            The qlen field is dependent on their being a self alignment in the alignment
            file, as then the qlen=tlen. If not present, qlen will be listed as 0.

            Note also the coverage is determined by alnlen/max(qlen, tlen)

            The output file is a .m8 file (e.g. tab delimited) and has the following
            columns:
            - query
            - target
            - query_id
            - target_id
            - alnlen
            - qlen
            - tlen
            - cov
            - pident
            - rmsd
            - z
            """
        ),
    )
    parser_aln_parse_dali.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        help="""
        Path to the DALI output file. The Summary output must be listed first!
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-s",
        "--structure_key",
        type=str,
        required=True,
        help="""
        Path to a tab-delimited file of format structure_name[delimiter]ID, where
        the ID is a 4-digit identifier used during the DALI alignment. This lets you
        convert the identifiers back. You can also import multiple separate structure
        key files! This entry should then be a comma-delimited list of paths to
        each structure key.
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output, now tab-delimited outfile. Will be written to (not appended)
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-z",
        "--min_z",
        type=float,
        required=False,
        default=0,
        help="""
        Alignments with a z-score below this number will be excluded. [Default: 0]
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-c",
        "--min_cov",
        type=float,
        required=False,
        default=0,
        help="""
        Alignments with a coverage below this number will be excluded. [Default: 0]
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-l",
        "--min_alnlen",
        type=float,
        required=False,
        default=0,
        help="""
        Alignments with an alnlen below this number will be excluded. [Default: 0]
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-r",
        "--min_rsmd",
        type=float,
        required=False,
        default=0,
        help="""
        Alignments with a rmsd below this number will be excluded. [Default: 0]
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-d",
        "--structure_key_delim",
        type=str,
        required=False,
        default=",,",
        help="""
        The tab delimiter of the structure_key. [Default: double comma]
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-n",
        "--query_name",
        type=str,
        required=False,
        default="",
        help="""
        A string indicating the name of the query, will go in the 'query' column of the
        output file. If this flag isn't used, the query will be searched for in the
        structure_key. If the query isn't in the structure_key either, the query column
        will be left blank.
        """,
    )
    parser_aln_parse_dali.add_argument(
        "-L",
        "--qlen_sheet",
        type=str,
        required=False,
        default="",
        help="""
        DALI doesn't report query len. Thus, if there is no self alignment, it is
        impossible to determine query len. This allows you to input a file of proteins
        and their associated length.
        
        Path to a file containing protein names and their associated amino acid sequence
        lengths. The first column should be the protein name, second column the length,
        and the file should be tab-delimited with one protein per line. If a protein's
        length is not present here, the script can determine qlen if there is a self
        alignment in the alignment file. Otherwise, qlen will be kept as 0.

        In this file, each file name can be the basename (with no .fasta or .pdb),
        or you can keep the .pdb or .fasta suffix.
        """,
    )
    parser_aln_parse_dali.set_defaults(func=call_aln_parse_dali_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_parse_dali subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_parse_dali_matrix = subparsers.add_parser(
        "aln_parse_dali_matrix",
        help=(
            """
            This subcommand takes in a DALI matrix file and/or a DALI dendogram files,
            and uses the specified key to convert each ID to its proper name.
            """
        ),
    )
    parser_aln_parse_dali_matrix.add_argument(
        "-k",
        "--key",
        type=str,
        required=True,
        help="""
        Path to a tab-delimited file of format structure_name,,ID, where
        the ID is a 4-digit identifier used during the DALI alignment. This lets you
        convert the identifiers back.
        """,
    )
    parser_aln_parse_dali_matrix.add_argument(
        "-t",
        "--tree",
        type=str,
        required=False,
        default="",
        help="""
        Path to the DALI tree file (ends in newick or newick_unrooted).
        """,
    )
    parser_aln_parse_dali_matrix.add_argument(
        "-T",
        "--tree_out",
        type=str,
        required=False,
        default="",
        help="""
        Path to the output tree file.
        """,
    )
    parser_aln_parse_dali_matrix.add_argument(
        "-m",
        "--matrix",
        type=str,
        required=False,
        default="",
        help="""
        Path to the DALI  similarity file.
        """,
    )
    parser_aln_parse_dali_matrix.add_argument(
        "-M",
        "--matrix_out",
        type=str,
        required=False,
        default="",
        help="""
        Path to the output, formatted matrix file.
        """,
    )
    parser_aln_parse_dali_matrix.set_defaults(func=call_aln_parse_dali_matrix_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_query_uniprot subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_query_uniprot = subparsers.add_parser(
        "aln_query_uniprot",
        help=(
            """
            This subcommand takes in a file that contains uniprotIDs (or alphafold IDs,
            even if they are formatted like e.g. AF-K0EZQ3-F1-model_v2.pdb.gz). You just
            need to specify the 0-indexed column position of the ID.

            This script uses the uniprot REST API to download information for each
            uniprotID. It offers the funcionality to save a cache containing the raw
            download information to prevent unncessary downloading.

            Besides making or updating a cache file, the main purpose of this script
            is to output a flat, tab-delimited file with the columns uniprotID,
            geneName, and fullName (fullName is a descriptive protein name). 
            """
        ),
    )
    parser_aln_query_uniprot.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="""
        Path to the input file that contains the uniprot IDs. If there are multiple
        columns, this must be tab-delimited.
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-o",
        "--uniprot_lookup_output",
        type=str,
        required=True,
        help="""
        This is the main output file. Will be a tab-delimited file with the columns
        uniprotID, geneName, fullName
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-c",
        "--infile_col",
        type=int,
        required=False,
        default=1,
        help="""
        Default: 1
        This is the 0-indexed column holding the uniprot of alphafold IDs.
        """,
    )
    parser_aln_query_uniprot.add_argument(
        "-u",
        "--uniprot_cache",
        type=str,
        required=False,
        default="",
        help="""
        A pkl file containing a cache from former uniprot downloads. This is optional.
        If specified, this file will be read in and will be updated with this script's
        downloads.
        """,
    )
    parser_aln_query_uniprot.set_defaults(func=call_aln_query_uniprot_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_taxa_counts subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_taxa_counts = subparsers.add_parser(
        "aln_taxa_counts",
        help=(
            """
            This takes in a cluster file (required columns are cluster_ID, cluster_rep,
            cluster_member, and cluster_count) and tallies up the taxons for each
            cluster. It makes a tidy file for each cluster where, for every taxon at
            every level, it specifies the count.

            The cluster file is assumed to be generated from an all-by-all alignment,
            perhaps with some additional merging steps. If you are also interested in
            adding taxonomy count information for the targets of a search of the
            cluster members against a separate database, you can enter an alignment
            file to this script. In the event an alignment file is provided, taxonIDs
            from the TARGET will be added to the cluster_ID of the QUERY.

            The output file has the
            following columns:
            cluster_ID, cluster_rep, cluster_count, superkingdom, level, taxon, count.
            """
        ),
    )
    parser_aln_taxa_counts.add_argument(
        "-c",
        "--cluster_file",
        type=str,
        required=True,
        help="""
        Path the cluster file. 
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path the output taxon counts file.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-F",
        "--cluster_file_fields",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited list of cluster file fields if none are present as headers in
        the cluster file. Only required if the cluster file doesn't come with headers.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=False,
        default="",
        help="""
        Path the alignment file.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-f",
        "--alignment_fields",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited list of alignment fields if none are present as headers in the
        alignment file.
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-t",
        "--taxonomy_levels",
        type=str,
        required=False,
        default="superkingdom,phylum,class,order,family,genus,species",
        help="""
        Taxonomy levels to count and output. Comma-delimited string of taxonomy levels.
        Default: superkingdom,phylum,class,order,family,genus,species
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-d",
        "--taxonID_finder_delimiter",
        type=str,
        required=False,
        default="__",
        help="""
        This script will try to find the taxonID in the query or target strings of each
        alignment (although if there is a taxid column in the alignment_object, it will
        use that taxonID for the target taxonID). To parse for the taxonID, this is
        the delimiter. [Default: __]
        """,
    )
    parser_aln_taxa_counts.add_argument(
        "-p",
        "--taxonID_finder_pos",
        type=int,
        required=False,
        default=-1,
        help="""
        This is where in the delimited string the taxonID is located. [Default: -1]
        """,
    )
    parser_aln_taxa_counts.set_defaults(func=call_parser_aln_taxa_counts_main)

    # -------------------------------------------------------------------------------- #
    # Parser for aln_trim_target_fasta subcommand
    # -------------------------------------------------------------------------------- #
    parser_aln_trim_target_fasta = subparsers.add_parser(
        "aln_trim_target_fasta",
        help=(
            """
        Trim the fasta sequences of target accessions from alignment files to the longest or shortest target alignment length.
        """
        ),
    )
    parser_aln_trim_target_fasta.add_argument(
        "-a",
        "--alignment_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to alignment file
        """,
    )
    parser_aln_trim_target_fasta.add_argument(
        "-c",
        "--alignment_fields",
        type=str,
        required=True,
        default="",
        help="""
        A string of column names for the alignment file
        """,
    )
    parser_aln_trim_target_fasta.add_argument(
        "-f",
        "--fasta_file",
        type=str,
        required=True,
        default="",
        help="""
        Path to fasta file
        """,
    )
    parser_aln_trim_target_fasta.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        default=1,
        help="""
        Name for output file
        """,
    )
    parser_aln_trim_target_fasta.add_argument(
        "-s",
        "--short_aln_type",
        type=int,
        required=False,
        default=1,
        help="""
        If 1, find the shortest alignment length.
        If 0, find the longest alignment length. 
        Default is 1. 
        """,
    )
    parser_aln_trim_target_fasta.set_defaults(func=call_aln_trim_target_fasta_main)

    # -------------------------------------------------------------------------------- #
    # Parser for tab_add_taxonomy subcommand
    # ----------------------------------------------------------------------------
    parser_tab_add_taxonomy = subparsers.add_parser(
        "tab_add_taxonomy",
        help=(
            """
            Adds taxonomy information to any file. There must be a column named
            'taxid', that contains each taxon ID. This file will essentially add
            additional columns, one per desired taxonomy level, with the corresponding
            taxon names.
            """
        ),
    )
    parser_tab_add_taxonomy.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="""
        Path to the input tabular file. Should be tab-delimited.
        """,
    )
    parser_tab_add_taxonomy.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="""
        Path to the output file. Will be tab-delimited.
        """,
    )
    parser_tab_add_taxonomy.add_argument(
        "-c",
        "--colnames",
        type=str,
        required=False,
        default="",
        help="""
        Comma-delimited string of the column names. If left blank, will use the first
        row of the infile as the colnames. For this script, one colname must be 'taxid'
        """,
    )
    parser_tab_add_taxonomy.add_argument(
        "-t",
        "--taxonomy_levels",
        type=str,
        required=False,
        default="superkingdom,phylum,class,order,family,genus,species",
        help="""
        Comma-delimited string taxonomy levels to be output.
        Default: 'superkingdom,phylum,class,order,family,genus,species'
        """,
    )
    parser_tab_add_taxonomy.set_defaults(func=call_parser_tab_add_taxonomy_main)

    # -------------------------------------------------------------------------------- #
    # Parser for plot_pae subcommand
    # -------------------------------------------------------------------------------- #
    parser_plot_pae = subparsers.add_parser(
        "plot_pae",
        help=(
            """
            This subcommand produces a PAE matrix plot when given a colabfold scores
            json file. The output file type is specified by the suffix of the out_image
            argument.
            """
        ),
    )
    parser_plot_pae.add_argument(
        "-s",
        "--scores",
        type=str,
        required=True,
        help="""
        Path to the colabfold PAE scores file, in json format.
        """,
    )
    parser_plot_pae.add_argument(
        "-o",
        "--out_image",
        type=str,
        required=True,
        help="""
        Path to the resultant PAE matrix image. The filetype of the output is determined
        by the suffix of this argument.
        """,
    )
    parser_plot_pae.set_defaults(func=call_plot_pae_main)

    # -------------------------------------------------------------------------------- #
    # Parser for seq_split_fasta subcommand
    # ----------------------------------------------------------------------------
    parser_seq_split_fasta = subparsers.add_parser(
        "seq_split_fasta",
        help=(
            """
            Tool to split a fasta file with multiple entries into individual fasta files 
            with one entry each.
            """
        ),
    )
    parser_seq_split_fasta.add_argument(
        "-i",
        "--in_fasta",
        type=str,
        required=True,
        help="""
        Path to the input fasta.
        """,
    )
    parser_seq_split_fasta.add_argument(
        "-o",
        "--outfile_dir",
        type=str,
        required=True,
        help="""
        Path to the output file directory
        """,
    )

    parser_seq_split_fasta.set_defaults(func=call_seq_split_fasta_main)

    # -------------------------------------------------------------------------------- #
    # Parser for seq_chunk subcommand
    # -------------------------------------------------------------------------------- #
    parser_seq_chunk = subparsers.add_parser(
        "seq_chunk",
        help=(
            """
            Tool to split sequences in an input fasta into overlapping or not
            overlapping chunks. Can write all resultant sequences to a single output
            file, or to separate output files (one per header). The header of
            chunks are >PART{N}_{header}.

            For example, if a sequence is 2300AA long but you desire sequences
            of 1000 max, this script can either generate the following:
            1) 1-1000, 1001-2000, 2001-2300. (if -v is NOT specified)
            2) 1-1000, 501-1500, 1001-2000, 1501-2300, 2001-2300 (if -v is specified)
            """
        ),
    )
    parser_seq_chunk.add_argument(
        "-i",
        "--in_fasta",
        type=str,
        required=True,
        help="""
        Path to the input fasta.
        """,
    )
    parser_seq_chunk.add_argument(
        "-o",
        "--out_fasta",
        type=str,
        required=True,
        help="""
        Path to the output fasta. If -n is specified, will put every fasta entry
        into a separate file and the entry in the -o switch will specify the
        base path - should be a directory in this case.
        """,
    )
    parser_seq_chunk.add_argument(
        "-m",
        "--max_seq_length",
        type=int,
        required=True,
        help="""
        The maximum sequence length.
        """,
    )
    parser_seq_chunk.add_argument(
        "-s",
        "--minimum_sequence_output_size",
        type=int,
        required=False,
        default=50,
        help="""
        Will not output any sequences below this size. This is because colabfold
        fails when sequences are very small. [50]
        """,
    )
    parser_seq_chunk.add_argument(
        "-v",
        "--overlapping_chunks",
        type=arg_str2bool,
        required=False,
        default=False,
        nargs="?",
        const=True,
        help="""
        If specified, each chunk will overlap by max_seq_length/2. If not
        specified, chunks will be just
        1:max_seq_length, max_seq_length+1:max_seq_length*2, etc
        """,
    )
    parser_seq_chunk.add_argument(
        "-n",
        "--individual",
        type=arg_str2bool,
        required=False,
        default=False,
        nargs="?",
        const=True,
        help="""
        If specified, each fasta entry will be passed to a separate file
        named by its header. Furthermore, it will assume that args.out_fasta is
        the base path.
        """,
    )
    parser_seq_chunk.set_defaults(func=call_seq_chunk_main)

    # -------------------------------------------------------------------------------- #
    # Parser for seq_multimerize subcommand
    # -------------------------------------------------------------------------------- #
    parser_seq_multimerize = subparsers.add_parser(
        "seq_multimerize",
        help=(
            """
            This subcommand combines input fastas to generate a multimierzed fasta
            containing :'s separating sequence. The cardinality of the input files can
            be specified to generate different kinds of homo- or hetero-complexes.
            Inputting the path to one file and specifing a cardinality of 2, for
            example, will generate an output file like the following:
            >args.header
            seq:seq
            """
        ),
    )
    parser_seq_multimerize.add_argument(
        "-i",
        "--infiles",
        type=str,
        required=True,
        help="""
        Path to one or more input fasta files. Each fasta file should have only one
        sequence. If there are multiple input files, this input is a *DOUBLE COMMA
        DELIMITED STRING*. E.g. an appropriate input for two files would be 
        '/path/to/file/file1.fasta,,/path/to/file/file2.fasta'
        """,
    )
    parser_seq_multimerize.add_argument(
        "-H",
        "--header",
        type=str,
        required=True,
        help="""
        This is the desired header of the output file.
        """,
    )
    parser_seq_multimerize.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="""
        Path to the output file.
        """,
    )
    parser_seq_multimerize.add_argument(
        "-c",
        "--cardinality",
        type=str,
        required=False,
        default="2",
        help="""
        This is a *DOUBLE COMMA DELIMITED STRING* of integers. There should be an equal
        number of intergers as there are paths in the infiles argument. Each integer
        value indicates the number of times the corresponding infile's sequence should
        be present in the output sequence.
        Example: '1,,2' will end up with the first sequence once and the second sequence
        twice.
        [Default: '2' (which makes a dimer given one input sequence)]
        """,
    )
    parser_seq_multimerize.set_defaults(func=call_seq_multimerize_main)

    # -------------------------------------------------------------------------------- #
    # Parser for seq_parse_genbank subcommand
    # -------------------------------------------------------------------------------- #
    parser_seq_parse_genbank = subparsers.add_parser(
        "seq_parse_genbank",
        help=(
            """
            This subcommand parses a genbank file (based on a nuclear accesion!) into
            an output fasta and an output table.

            The output fasta will have headers with the following information:
            {genome_acc}{args.delimiter}{protein_id}{args.delimiter}{locus_tag}{args.delimiter}{protein_order}
            This is equivelant to the "output_name"

            The output table is CSV FORMATTED, with the following columns:
            {output_name},{genome_acc},{locus_tag},{protein_id},{start},{end},{strand},{protein_order},{organism_name},{protein_name}

            Note that if you desire to only process a subset of genbank entires, you
            can provide a file with the genome accessions (no version!) that you
            desire.
            """
        ),
    )
    parser_seq_parse_genbank.add_argument(
        "-i",
        "--gb",
        type=str,
        required=True,
        help="""
        Path to the input genbank file. This should be a nucleotide genbank.
        """,
    )
    parser_seq_parse_genbank.add_argument(
        "-f",
        "--out_fasta",
        type=str,
        required=True,
        help="""
        Path to the output fasta.
        """,
    )
    parser_seq_parse_genbank.add_argument(
        "-c",
        "--out_csv",
        type=str,
        required=True,
        help="""
        Path to the output csv.
        """,
    )
    parser_seq_parse_genbank.add_argument(
        "-d",
        "--delimiter",
        type=str,
        required=False,
        default="__",
        help="""
        This is the delimiter that separates all fields of the fasta header. Default
        is double underscore "__"
        """,
    )
    parser_seq_parse_genbank.add_argument(
        "-n",
        "--nuc_accs_to_keep",
        type=str,
        required=False,
        default="",
        help="""
        Path to a file that contains desired nucleotide accessions, one per line. 
        Only genbank entries with that accession will be parsed. Note that the accession
        should not have the version (e.g. YP_324124 works, YP_324124.1 does not work).

        This argument is optional. If not specified, all genbank entires will be parsed.
        """,
    )

    parser_seq_parse_genbank.set_defaults(func=call_parser_seq_parse_genbank)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_detect_interaction subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_detect_interaction = subparsers.add_parser(
        "struc_detect_interaction",
        help=(
            """
            This subcommand takes in a structure predction (PDB file) and its associated
            PAE file (from colabfold) that was generated with AF Multimer between
            two molecules. Thus, the structure prediction should have two chains.

            This script does two main things:
            1) Determins which residues of chain 1 are in contact with residues from
            chain 2. Residues are considered in contact if the distance between an atom
            in a chain 1 residue and an atom in a chain 2 residue is less than the sum
            of their van der Waals radii plus 0.5 angstroms. This script counts both the total
            number of pairwise interactions as well as the total number of residues
            (from either chain) that are at the interface of the two chains.
            2) Using the PAE matrix, clusters the residues and determins if there is a
            cluster that contains residues from both chains.

            The output file is tab-delimited with the following columns:
            - member1
            - member2
            - average PAE of residues that are interacting across chains
            - number of interactions between residues in each chain
            - total number of residues at the interface of the two chains
            - number of residues in chain1
            - number of residues in chain2
            """
        ),
    )
    parser_struc_detect_interaction.add_argument(
        "-s",
        "--structure",
        type=str,
        required=True,
        help="""
        Path to the structure file
        """,
    )
    parser_struc_detect_interaction.add_argument(
        "-p",
        "--pae",
        type=str,
        required=True,
        help="""
        Path to the structure file
        """,
    )
    parser_struc_detect_interaction.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="""
        Path to the output file, which is tab delimited and has the following columns:
        - member1
        - member2
        - interaction (True or False)
        - # of residues in chain1 that have a C-alpha within distance_cutoff
            angstroms of a C-alpha from chain2.
        - # of residues in chain2 that have a C-alpha within distance_cutoff
            angstroms of a C-alpha from chain1.
        """,
    )
    parser_struc_detect_interaction.add_argument(
        "-D",
        "--delimiter",
        type=str,
        required=False,
        default="__",
        help="""
        Delimiter present in the structure file basename that separates the molecules
        that were folded together. For example, a structure file of
        /path/to/file/mol1__mol2.pdb with a delimiter of __ will identify mol1 and mol2
        as the member names.
        [Default: '__']
        """,
    )
    parser_struc_detect_interaction.set_defaults(
        func=call_struc_detect_interaction_main
    )

    # -------------------------------------------------------------------------------- #
    # Parser for struc_disorder subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_disorder = subparsers.add_parser(
        "struc_disorder",
        help=(
            """
            This takes an input structure and calculates the number of residues that
            are considered ordered, disordered, or intermediate. A residue is considered
            ordered if it is in a stretch of at least n_sequential residues that have a
            pLDDT of >= order_cutoff. A residue is considered disordered if it is in a
            stretech of at least n_sequential residues <= disorder_cutoff.

            This returns an output file with the following columns:
            - basename of the input structure
            - number of ordered residues
            - number of disordered residues
            - number of intermediate residues (neither ordered or disordered)
            - total number of residues
            - there_is_a_domain: yes or no. This checks that there is at least one
            stretech of continuous residues that have ordered pLDDTs. The required
            stretch size is args.check_for_domain_len
            """
        ),
    )
    parser_struc_disorder.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file
        """,
    )
    parser_struc_disorder.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="""
        Path to the output file.
        """,
    )
    parser_struc_disorder.add_argument(
        "-c",
        "--disorder_cutoff",
        type=int,
        required=False,
        default=50,
        help="""
        Integer value. Residues that have a pLDDT <= this value and are in a stretch of
        n_sequential residues are considered disordered. [Default: 50]
        """,
    )
    parser_struc_disorder.add_argument(
        "-C",
        "--order_cutoff",
        type=int,
        required=False,
        default=60,
        help="""
        Integer value. Residues that have a pLDDT >= this value and are in a stretch of
        n_sequential residues are considered ordered. [Default: 60]
        """,
    )
    parser_struc_disorder.add_argument(
        "-n",
        "--n_sequential",
        type=int,
        required=False,
        default=6,
        help="""
        Integer value. There must be n_sequential residues for a residue to be called
        ordered or disordered. [Default: 6]
        """,
    )
    parser_struc_disorder.add_argument(
        "-d",
        "--check_for_domain_len",
        type=int,
        required=False,
        default=50,
        help="""
        Integer value. This checks that there is at least one well-folded stretch of
        residues with at least this pLDDT. [Default: 50]
        """,
    )
    parser_struc_disorder.set_defaults(func=call_struc_disorder)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_download subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_download = subparsers.add_parser(
        "struc_download",
        help=(
            """
            This subcommand takes in a file of uniprot IDs and downloads the
            AF2 database pdb and pae files to the indicated directory. Furthermore,
            if any additional information is present in the tabular infile it will be
            appended to the output files - this is a good way to lable the files with
            information like taxonomyID, etc.
            """
        ),
    )
    parser_struc_download.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="""
        Path to a file containing uniprotIDs. This file can just have the uniprotIDs
        as a single column, OR it can be tab-delimited and have additional information
        for each uniprotID on subsequent columns (e.g. taxonomy ID, etc). If
        this file has multiple columns, you must set the infile_columns parameter
        accordingly. If there are multiple columns, the information from those columns
        will be appended to the end of the pdb and pae file outputs via a delimted
        field specified by additional_field_delimiter.
        """,
    )
    parser_struc_download.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="""
        Path to the output directory that will contain all output pdb and pae files.
        """,
    )
    parser_struc_download.add_argument(
        "-c",
        "--infile_columns",
        type=str,
        required=False,
        default="uniprotID",
        help="""
        Default: uniprotID.\n
        This is a comma-delimeted string that indicates the names of the columns present
        in the input file. If there is only one column of uniprotIDs you can leave it as
        default.
        """,
    )
    parser_struc_download.add_argument(
        "-d",
        "--additional_field_delimiter",
        type=str,
        required=False,
        default="__",
        help="""
        This is the delimiter that will be used when adding the additional fields if
        present in the infile.
        """,
    )
    parser_struc_download.set_defaults(func=call_struc_download_main)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_extract_chains subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_extract_chains = subparsers.add_parser(
        "struc_extract_chains",
        help=(
            """
            This subcommand extracts one or more chains from an input structure,
            and writes them to a new pdb file. The desired chains should be input as a
            comma-delimited string.
            """
        ),
    )
    parser_struc_extract_chains.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the input structure file.
        """,
    )
    parser_struc_extract_chains.add_argument(
        "-c",
        "--chains",
        type=str,
        required=True,
        help="""
        Comma-delimited list containing the chainID(s) of the chains you want to write
        to the output file
        """,
    )
    parser_struc_extract_chains.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="""
        Path to the output structure file containing only the specified chain(s).
        """,
    )
    parser_struc_extract_chains.set_defaults(func=call_struc_struc_extract_chains_main)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_find_motif subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_find_motif = subparsers.add_parser(
        "struc_find_motif",
        help=(
            """
            Given a motif of structure [OPTIONS]xxx[OPTIONS]xx, where x indicates
            any amino acid and [] indicate any of the amino acids present within the
            brackets, this returns the match and position start/end of the motif
            present in the input sequence.

            The input can be a structure file, a fasta, or just a sequence. The output
            is tab-delimited and printed to the screen, with the columns
            - match
            - start
            - end

            Where start and end are 1-indexed.
            """
        ),
    )
    parser_struc_find_motif.add_argument(
        "-m",
        "--motif",
        type=str,
        required=True,
        help="""
        String indicating the motif you're searching for. If an AA at a position doesn't
        matter, use a lowercase x. If it does matter, put an upercase AA single letter
        digit. If a position can have one of multiple possible AAs, use brackets... e.g.
        [RK] indicates a position can either be R or K.

        Make sure to wrap this in quotes when inputting.
        """,
    )
    parser_struc_find_motif.add_argument(
        "-s",
        "--structure",
        type=str,
        required=False,
        default="",
        help="""
        Path to the structure file.
        """,
    )
    parser_struc_find_motif.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=False,
        default="",
        help="""
        Path to the fasta file.
        """,
    )
    parser_struc_find_motif.add_argument(
        "-q",
        "--seq",
        type=str,
        required=False,
        default="",
        help="""
        A string with the AA sequence.
        """,
    )
    parser_struc_find_motif.set_defaults(func=call_struc_find_motif)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_get_contact_probability subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_get_contact_probability = subparsers.add_parser(
        "struc_get_contact_probability",
        help="""
        This script calculates the internal contact probability between two chains
        folded using colabfold. This can be from any model type (e.g. 
        alphafold2_multimer_v3, or even alphafold2_ptm) as long as two chains were
        predicted using colabfold. 

        This script parses the entire colabfold output directory. For a given model
        type, if different models were run (between 1-5) this script will parse 
        every model. 

        This will measure the highest contact probability in the submatrix of residues
        that cross chains A and B. The user can set a distance cutoff - if set to e.g.
        8, this script, for every residue-residue pair between chains, will determine
        the probability the residue pairs are within 8 angstroms. The highest 
        residue-residue contact probability is the result. I've annotated the functions
        in this script if you want to learn more about the technical details.

        For each model, there are required input files:
        1) The a3m file. This is only present once regardless of the number of models
            run. The critical part of this file is the first line, which lists the
            length of the two protein chains.
        2) The pickle file. This is the ouput pickle object from each model. This is
            the source of the logits from which we calculate interaction probability.
        3) The json file. This is the scores json file that lists PAE, but also at the
            end contains the iptm.

        The output file is tab delimited and has the following columns, with one 
        line per model:
        - sample_name (the prefix of the input files)
        - model type (e.g. alphafold2_ptm)
        - rank (1-5, ranked based on however colabfold was ranking structures. Usually
            pLDDT. E.g. rank 1 is the structure with highest pLDDT)
        - model (1-5, this is the model run)
        - iptm (interaction pTM score reported by colabfold)
        - contact probability (highest cross-chain residue-residue contact probability
            as discussed above)
        """,
    )
    parser_struc_get_contact_probability.add_argument(
        "-i",
        "--in_dir",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input directory. This should be produced by colabfold, and must 
        have one or more .pickle files (output using setting --save-all) in addition
        to one a3m file and one or more json files.
        """,
    )
    parser_struc_get_contact_probability.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        default="",
        help="""
        Output file, containing one line per model present in the input directory.
        Columns are as follows:
        - sample_name (the prefix of the input files)
        - model type (e.g. alphafold2_ptm)
        - rank (1-5, ranked based on however colabfold was ranking structures. Usually
            pLDDT. E.g. rank 1 is the structure with highest pLDDT)
        - model (1-5, this is the model run)
        - iptm (interaction pTM score reported by colabfold)
        - contact probability (highest cross-chain residue-residue contact probability
            as discussed)
        """,
    )
    parser_struc_get_contact_probability.add_argument(
        "-d",
        "--distance_cutoff",
        type=int,
        required=False,
        default=8,
        help="""
        [Default: 8] Number of angstroms under which we are extracting the contact
        probability. See the module notes for more information.
        """,
    )

    parser_struc_get_contact_probability.set_defaults(
        func=call_struc_get_contact_probability
    )

    # -------------------------------------------------------------------------------- #
    # Parser for struc_get_domains subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_get_domains = subparsers.add_parser(
        "struc_get_domains",
        help="""
        Extract domains from a structure using chainsaw domain predictions. Notably, this script is
        designed for pdb files. 
        """,
    )
    parser_struc_get_domains.add_argument(
        "-s",
        "--structure_file_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input structure in .pdb format.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-c",
        "--chainsaw_file_path",
        type=str,
        required=True,
        default="",
        help="""
        Path to the chainsaw.txt file. The first row contains column names.
        The chainsaw file can consist of a single chainsaw output or be a concatenation of multiple chainsaw outputs.
        """,
    )

    parser_struc_get_domains.add_argument(
        "-m",
        "--min_domain_length",
        type=int,
        required=False,
        default=1,
        help="""
        Optional minimum length of each domain. Default value is 1. Domains that do not meet the minimum length requirment
        will be excluded and not output a pdb file.
        """,
    )
    parser_struc_get_domains.add_argument(
        "-o",
        "--outfile_dir",
        type=str,
        required=True,
        default="",
        help="""
        Directory of the output domains. Files will be labeled
        {output_dir}/{basename of structure}_domain_{domain_start}_{domain_end}.pdb.
        Note that the domain number will be 1-indexed.
        """,
    )

    parser_struc_get_domains.set_defaults(func=call_struc_get_domains)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_get_itpm subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_get_iptm = subparsers.add_parser(
        "struc_get_iptm",
        help="""
        Extract the iPTM value from the colabfold json. Appends to an output file with
        format [json_basename]\t[iptm]\n. 
        """,
    )
    parser_struc_get_iptm.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        default="",
        help="""
        Path to the input json file from colabfold. 
        """,
    )
    parser_struc_get_iptm.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        default="",
        help="""
        Path to the output tsv - will have the following columns:
        1. basename of the input json
        2. iptm value

        This output file will be appended to if it already exists
        """,
    )

    parser_struc_get_iptm.set_defaults(func=call_struc_get_iptm)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_qc subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_qc = subparsers.add_parser(
        "struc_qc",
        help=(
            """
            Given a structure, determines the percentage of residues that have at least
            the specified pLDDT. The output is returned to STDOUT!! It is tab-delimited
            and has the following columns:
            - structure_file (the basename of the file, including suffix)
            - number of residues
            - number of residues that pass the pLDDT threshold
            - proportion of residues that pass the pLDDT threshold (this will be a
                decimal between 0 and 1)
            """
        ),
    )
    parser_struc_qc.add_argument(
        "-s",
        "--structure_file_path",
        type=str,
        required=True,
        help="""
        Path to the structure file
        """,
    )
    parser_struc_qc.add_argument(
        "-p",
        "--plddt_cutoff",
        type=int,
        required=True,
        help="""
        The minimum pLDDT for a residue to be considered passing. This is an integer
        value between 0 and 100.
        """,
    )
    parser_struc_qc.set_defaults(func=call_parser_struc_qc)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_rebase subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_rebase = subparsers.add_parser(
        "struc_rebase",
        help=(
            """
            Simple subcommand that renumbers all residues in a structure such that
            the first residue is #1 and all residues are sequential (e.g. it takes out
            numeric gaps in residue numbers).
            """
        ),
    )
    parser_struc_rebase.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file in pdb format.
        """,
    )
    parser_struc_rebase.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=True,
        help="""
        Path to the output structure file.
        """,
    )
    parser_struc_rebase.set_defaults(func=call_struc_rebase_main)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_to_plddt subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_to_plddt = subparsers.add_parser(
        "struc_to_plddt",
        help=(
            """
            Simple subcommand that returns the average plddt of the input structure
            file. If --out_file is not specified, the average plddt is simply printed
            to the screen. If --out_file is specified, the output file will be
            APPENDED to with the following:
            [basename input structure_file]\\t[plddt]\\n
            """
        ),
    )
    parser_struc_to_plddt.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file in pdb format.
        """,
    )
    parser_struc_to_plddt.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=False,
        default="",
        help="""
        Path to a file the sequence will be appended to. If left blank, sequence will be
        printed to the screen.
        """,
    )
    parser_struc_to_plddt.set_defaults(func=call_struc_to_plddt)

    # -------------------------------------------------------------------------------- #
    # Parser for struc_to_seq subcommand
    # -------------------------------------------------------------------------------- #
    parser_struc_to_seq = subparsers.add_parser(
        "struc_to_seq",
        help=(
            """
            Simple subcommand that returns the amino-acid sequence of a specified
            structure (in pdb format). The AA sequence will be APPENDED to the outfile
            if specified, or printed to the screen if -o --out_file is not specified.
            """
        ),
    )
    parser_struc_to_seq.add_argument(
        "-s",
        "--structure_file",
        type=str,
        required=True,
        help="""
        Path to the structure file in pdb format.
        """,
    )
    parser_struc_to_seq.add_argument(
        "-o",
        "--out_file",
        type=str,
        required=False,
        default="",
        help="""
        Path to a file the sequence will be appended to. If left blank, sequence will be
        printed to the screen.
        """,
    )
    parser_struc_to_seq.add_argument(
        "-H",
        "--header",
        type=str,
        required=False,
        default="",
        help="""
        Header of the entry if writing to a fasta. Only required if -o --out_file is
        specified.
        """,
    )
    parser_struc_to_seq.set_defaults(func=call_struc_to_seq)

    # ----------------------------------------------------------------------------------#
    # Parse the args and call the function associated with the subcommand
    # ----------------------------------------------------------------------------------#
    args = parser.parse_args()
    args.func(args)


def call_aln_add_clusters(args):
    from scripts.aln_add_clusters import aln_add_clusters_main

    aln_add_clusters_main(args)


def call_aln_add_uniprot_main(args):
    from scripts.aln_add_uniprot import (
        aln_add_uniprot_main,
    )

    aln_add_uniprot_main(args)


def call_aln_connection_map(args):
    from scripts.aln_connection_map import aln_connection_map_main

    aln_connection_map_main(args)


def call_aln_cluster_connected_component_main(args):
    from scripts.aln_cluster_connected_component import aln_cluster_connected_component_main

    aln_cluster_connected_component_main(args)


def call_aln_dali_alignment_attributes_main(args):
    from scripts.aln_dali_alignment_attributes import aln_dali_alignment_attributes_main

    aln_dali_alignment_attributes_main(args)


def call_aln_dali_motif_finder_main(args):
    from scripts.aln_dali_motif_finder import aln_dali_motif_finder_main

    aln_dali_motif_finder_main(args)


def call_aln_dali_to_pairwise_fastas_main(args):
    from scripts.aln_dali_to_pairwise_fastas import aln_dali_to_pairwise_fastas_main

    aln_dali_to_pairwise_fastas_main(args)


def call_aln_filter_main(args):
    from scripts.aln_filter import aln_filter_main

    aln_filter_main(args)


def call_aln_cigar_to_cov_main(args):
    from scripts.aln_cigar_to_cov import aln_cigar_to_cov_main

    aln_cigar_to_cov_main(args)


def call_aln_cluster_greedy_main(args):
    from scripts.aln_cluster_greedy import aln_cluster_greedy_main

    aln_cluster_greedy_main(args)


def call_aln_merge_clusters_main(args):
    from scripts.aln_merge_clusters import aln_merge_clusters_main

    aln_merge_clusters_main(args)


def call_aln_merge_main(args):
    from scripts.aln_merge import aln_merge_main

    aln_merge_main(args)


def call_aln_parse_dali_main(args):
    from scripts.aln_parse_dali import aln_parse_dali_main

    aln_parse_dali_main(args)


def call_aln_parse_dali_matrix_main(args):
    from scripts.aln_parse_dali_matrix import aln_parse_dali_matrix_main

    aln_parse_dali_matrix_main(args)


def call_aln_query_uniprot_main(args):
    from scripts.aln_query_uniprot import aln_query_uniprot_main

    aln_query_uniprot_main(args)


def call_aln_trim_target_fasta_main(args):
    from scripts.aln_trim_target_fasta import aln_trim_target_fasta_main

    aln_trim_target_fasta_main(args)


def call_parser_aln_taxa_counts_main(args):
    from scripts.aln_taxa_counts import aln_taxa_counts_main

    aln_taxa_counts_main(args)


def call_parser_seq_parse_genbank(args):
    from scripts.seq_parse_genbank import seq_parse_genbank_main

    seq_parse_genbank_main(args)


def call_parser_struc_qc(args):
    from scripts.struc_qc import struc_qc_main

    struc_qc_main(args)


def call_parser_tab_add_taxonomy_main(args):
    from scripts.tab_add_taxonomy import tab_add_taxonomy_main

    tab_add_taxonomy_main(args)


def call_plot_pae_main(args):
    from scripts.plot_pae import plot_pae_main

    plot_pae_main(args)


def call_seq_chunk_main(args):
    from scripts.seq_chunk import seq_chunk_main

    seq_chunk_main(args)


def call_seq_multimerize_main(args):
    from scripts.seq_multimerize import seq_multimerize_main

    seq_multimerize_main(args)


def call_seq_split_fasta_main(args):
    from scripts.seq_split_fasta import seq_split_fasta_main

    seq_split_fasta_main(args)


def call_struc_detect_interaction_main(args):
    from scripts.struc_detect_interaction import struc_detect_interaction_main

    struc_detect_interaction_main(args)


def call_struc_disorder(args):
    from scripts.struc_disorder import struc_disorder_main

    struc_disorder_main(args)


def call_struc_download_main(args):
    from scripts.struc_download import struc_download_main

    struc_download_main(args)


def call_struc_find_motif(args):
    from scripts.struc_find_motif import struc_find_motif_main

    struc_find_motif_main(args)


def call_struc_get_contact_probability(args):
    from scripts.struc_get_contact_probability import struc_get_contact_probability_main

    struc_get_contact_probability_main(args)


def call_struc_get_domains(args):
    from scripts.struc_get_domains import struc_get_domains_main

    struc_get_domains_main(args)


def call_struc_get_iptm(args):
    from scripts.struc_get_iptm import struc_get_iptm_main

    struc_get_iptm_main(args)


def call_struc_rebase_main(args):
    from scripts.struc_rebase import struc_rebase_main

    struc_rebase_main(args)


def call_struc_struc_extract_chains_main(args):
    from scripts.struc_extract_chains import struc_extract_chains_main

    struc_extract_chains_main(args)


def call_struc_to_plddt(args):
    from scripts.struc_to_plddt import struc_to_plddt_main

    struc_to_plddt_main(args)


def call_struc_to_seq(args):
    from scripts.struc_to_seq import struc_to_seq_main

    struc_to_seq_main(args)


# Keep these buffer lines here
#
#
if __name__ == "__main__":
    main()
