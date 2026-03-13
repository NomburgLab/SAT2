# SAT
Structural Analysis Toolkit - A python package for manipulation of structural data and structural alignments.  

![Tests](https://github.com/jnoms/SAT2/actions/workflows/main.yml/badge.svg)

# Installation
There are two methods to install this package - in a poetry environment, or within a conda environment.

## Installation in a poetry environment
1. Make sure the [Poetry package manager](https://python-poetry.org/) is installed.
1. Clone this repository
2. Enter the SAT directory, and install the package via poetry. This will download all dependencies into a poetry virtual environment.
```
poetry install
```
3. Now, you can run SAT in the following way:
```
poetry run sat.py <subcommand>
```

## Installation in a conda environment
1. Make sure the [conda package manager](https://docs.conda.io/en/latest/miniconda.html) is installed.
1. Clone this repository
2. Create a conda environment that contains the poetry package manager
```
conda create --name SAT -c conda-forge poetry=1.8.3 python=3.10.5
```
3. Activate environment. Enter the SAT directory and download the dependencies using poetry. Dependencies will be downloaded specifically into that conda environment.
```
conda activate SAT
poetry install
```
4. The SAT conda environment will now contain all dependencies. You can run SAT in the following way:
```
conda activate SAT #if conda enviornment is not active  

sat.py <subcommand>
```

## Test installation
Navigate to the sat directory and enter `pytest` (if using a conda environment) or `poetry run pytest` if using a poetry environment. This will make sure that all tests pass and sat is properly installed.  

## Note on ETE3 
When you run the tests or the first time you run any taxonomy-related script, ete3 will download a taxonomy database to ~/.etetoolkit/. **This database is 560MB** as of October 2022. Capability to specify the database download location is on the to-do list, but in the interim you can make a symlink from ~/.etetoolkit/ to wherever you want the database to reside (see below). This is particularly important if taxonomy related queries are going slowly, as that probably means your home directory has slow IO so you should symlink to somewhere with faster IO.    
```ln -s /desired/ete/database/location ~/.etetoolkit```  


# List of subcommands

## Structure-focused
`sat.py struc_clean` - Cleans a pdb file from water molecules, ions and other ligands, leaving only the proteins. 
`sat.py struc_detect_interaction` - For a co-folded prediction of two molecules, determine if the PAE matrix clusters across the molecules and suggests a potential interaction.  
`sat.py struc_disorder` - Get information on the number of residues in an input structure that are considered disordered and ordered.  
`sat.py struc_extract_chains` - Given an input structure file with multiple chains, write a new file with only the specified chain(s).  
`sat.py struc_find_motif` - Checks if there is a motif in a structure or sequence input.  
`sat.py struc_get_contact_probability` - Determines the probability that two proteins are interacting. This will probably replace struc_detect_interaction.  
`sat.py struc_get_domains` - Uses chainsaw predicted domain boundaries to extract domains from structure files.  
`sat.py struc_get_iptm` - Extract the iPTM value from a colanbfold json file.    
`sat.py struc_qc` - Get information on the fraction of residues that are at least a specified pLDDT - this can be good for filtration.  
`sat.py struc_rebase` - Rebases an input structure such that the first residue is residue #1, and all subsequent residues are sequential (e.g. removes numeric gaps present in discontinuous domains).  
`sat.py struc_to_plddt` - Prints the average pLDDT of a structure to the screen or appends to a specified file.  
`sat.py struc_to_seq` - Prints the primary amino acid sequence of a structure to the screen or appends to a specified file in fasta format.  


## Alignment-focused
`sat.py aln_add_clusters` - Adds foldseek clustering information to the foldseek tabular alignment file.  
`sat.py aln_add_uniprot` - After retreiving the uniprot unformation using aln_query_uniprot, adds the information as columns to the alignment file.  
`sat.py aln_cigar_to_cov` - Adds CIGAR-derived coverage columns (cigar_qcov, cigar_tcov) to an alignment file by parsing the CIGAR string.  
`sat.py aln_cluster_connected_component` - Generates connected component clusters from an alignment file. All query-target pairs that are connected will be placed into the same cluster.  
`sat.py aln_cluster_greedy` - Performs greedy set cover clustering, similar to foldseek/mmseqs cluster mode 0. Does not transitively connect members.  
`sat.py aln_connection_map` - This takes a cluster file (that has taxonomy information) and reports, for every pair of families, the number of clusters that they share.  
`sat.py aln_dali_alignment_attributes` -  Generates a csv file of aligned targets and queries and their attributes, such as qstart, qend, tstart, tend, etc.  
`sat.py aln_dali_motif_finder` - This parses the 'alignments' block of a DALI output file, and checks if a given residue or motif is present in each target at an indicated position of the structural alignment.  
`sat.py aln_dali_to_pariwise_fastas` - This subcommand converts the alignment fields of a DALI output file to pairwise fasta files, one per alignment. The usecase here is typically when you are trying to generate a structure-guided MSA.  
`sat.py aln_filter` - This filters for alignments below/above a specified value in a specified column, and can also filter to keep a maximum number of queries per alignment.  
`sat.py aln_merge` - This merges two alignment files.  
`sat.py aln_merge_clusters` - Merges two cluster files by adding a higher-level clustering to a nested cluster file.  
`sat.py aln_parse_dali` - This parses a Dalilite alignment output into a tab-delimited format. It can also filter based on various alignment statistics.  
`sat.py aln_parse_dali_matrix` - This parses a Dalilite matrix output - basically uses the DALI key to annotate it.  
`sat.py aln_query_uniprot` - Lets you look up alphafold or uniprot IDs using the Uniprot REST API, and get the geneName and fullName (an informative protein name) for each.  
`sat.py aln_taxa_counts` - Returns counts at desired taxonomic levels within each foldseek cluster.  
`sat.py aln_trim_target` - Trim the fasta sequences of target accessions from alignment files to the longest or shortest target alignment length.

## Sequence-focused  
`sat.py seq_chunk` - Splits a fasta file into overlapping or non-overlapping chunks.  
`sat.py seq_multimerize` - Combines one or more fasta sequences, separated by :'s, to be used for multimer prediction. Cardinality can be specified, so this is good to make any number of homo- and hetero-complexes.  
`sat.py seq_parse_genbank` - Parses a nucleotide genbank file into a fasta of proteins, as well as a convenient table.  
`sat.py seq_split_fasta` - Parses a fasta file with multiple entries into multiple fasta files with one entry each.      

## Plotting-focused  
`sat.py plot_pae` - Plots the colabfold PAE scores json file.    

## Manipulation of tabular files  
`sat.py tab_add_taxonomy` - If there is a column in a file with taxonomy IDs, this script looks up the lineage and adds the lineage to each line as a column.   

# SAT aln_add_uniprot
This script adds the uniprot information garther from aln_query_uniprot to a foldseek alignment file.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_add_uniprot -h`](.github/img/aln_add_uniprot.png)  

# SAT aln_cigar_to_cov
This subcommand adds CIGAR-derived coverage columns to an alignment file. It parses the CIGAR string to count the number of M (match) operations, then calculates:
- `cigar_qcov` = #Ms / qlen (query coverage)
- `cigar_tcov` = #Ms / tlen (target coverage)

These columns are inserted directly after the `cigar` column in the output.

The input alignment file must have `qlen`, `tlen`, and `cigar` columns. Column names can be auto-detected if the first line starts with "query", or provided via `--colnames` as a comma-delimited list.

**Example:**

Input:
```
query   target  qlen    tlen    cigar   evalue
q1      t1      100     200     80M     1e-10
```

Output:
```
query   target  qlen    tlen    cigar   cigar_qcov  cigar_tcov  evalue
q1      t1      100     200     80M     0.800       0.400       1e-10
```

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_cigar_to_cov -h`](.github/img/aln_cigar_to_cov.png)  

# SAT aln_cluster_connected_component
This subcommand generates connected component clusters from an alignment file. All query-target pairs that are connected (directly or transitively) will be placed into the same cluster. This is similar to foldseek/mmseqs cluster mode 1 (connected-component clustering).

The output is a foldseek-style cluster file with two columns: `cluster_rep` and `cluster_member`. The cluster representative is selected randomly from each cluster.

Column names can be auto-detected if the first line of the alignment file starts with "query". Otherwise, provide column names via `--colnames` as a comma-delimited list (must include "query" and "target").

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_cluster_connected_component -h`](.github/img/aln_cluster_connected_component.png)  

# SAT aln_cluster_greedy
This subcommand performs greedy set cover clustering, similar to foldseek/mmseqs cluster mode 0.

**Algorithm:**
1. Sort members by number of alignments (descending)
2. Pick the top member as cluster rep, assign all its alignments to its cluster
3. Remove assigned members from consideration
4. Repeat until all members are clustered
5. Perform a reassignment step where members are moved to a different cluster if they have a better alignment score to that cluster's representative

**Key difference from connected-component clustering:** This method does NOT transitively connect members. If A aligns to B and B aligns to C, but A doesn't align to C, then A and C may end up in different clusters. This typically produces more, smaller clusters than connected-component clustering.

The output is a foldseek-style cluster file with two columns: `cluster_rep` and `cluster_member`.

Use `--no_reassign` to skip the reassignment step. Use `-s/--score_column` to specify which column contains alignment scores (default: `alntmscore`).

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_cluster_greedy -h`](.github/img/aln_cluster_greedy.png)  

# SAT aln_connection_map
This subcommand takes in a cluster file that has taxonomy information (critically - family) and determines, for each pair of families, how many clusters exist in which both families have a member.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_aln_connection_mapecod_purity -h`](.github/img/aln_connection_map.png)  

# SAT aln_dali_alignment_attributes
This subcommands takes in a DALI alignment file (which must have an alignments field) and a key to generate a csv file of aligned targets and queries and their attributes, such as qstart, qend, tstart, tend, etc.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_dali_alignment_attributes -h`](.github/img/aln_dali_alignment_attributes.png)  

# SAT aln_dali_motif_finder
This subcommand parses a DALI alignments field (specified by the 'alignments' output format) and determines, for each target, if a specified motif or series of motifs is present at specified locations. The output of this script is a a file with two columns - target and target_id - for those targets that contain all indicated motifs. Note that adding a DALI key is optional, in which case 'target' will be blank. See below for detailed instructions about how to input motifs.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_dali_motif_finder -h`](.github/img/aln_dali_motif_finder.png) 


## Description of the motif_list input argument
*Important Note*: The motif query position (e.g. POS, see below) is based on the residues actually included in the structure. E.g. if you are looking for K26, but the PDB file started at residue 2, you need to use 25 as the value for POS.  

This is a complicated string with many sublists. I will describe them iterateively.  
        motif_list: MOTIF+MOTIF  
        MOTIF: POS_RESIDUES_FLEXIBILITY  
        - POS: a 1-indexed position of a residue in the query. This script will determine where in the alignment that position occurs, and then check for 
            the indicated motif (specified in RESIDUES) starting at that alignment
            position.  
        - RESIDUES: a comma-delimited list of residues. E.g., the length of this list 
            is equal to the length of the desired motif. For each individual motif
            position, multiple residue options are specified with forward-slash-
            delimited list. An 'X' indicates that E.g. R/H,X,H is searching for a motif
            that starts with R or H, followed by any residue, followed by H.  
        - FLEXIBILITY: an integer value indicating how flexible the motif positioning
            is. If set to 0, this means that the RESIDUES must start exactly at the
            alignment position derived from POS. If 1, the motif can start 1 residue
            beforehand or 1 residue after that position. And so on.   

        Examples of correct inputs, and their interpretations:  
        16_K_0                  K at the alignment position of query residue 16.  
        16_K_0+82_K_0           K at the alignment position of query residue 16 and 82.  
        16_K_1                  K in the target either at the alignment position 
                                corresponding to query reisdue 16, or one residue 
                                beforehand, or one residue afterwards.  
        72_R/H,X,K_1+103_H_0    Two motifs. First is looking for an R or H, followed by
                                any residue, followed by K. The position of this motif
                                should be at the alignment position of query reisdue 72
                                or one residue before or after that position. The second
                                motif is looking for an H at the alignment position 
                                corresponding to query residue 103.  

        Description of the POS (indicated by the input motif) vs the alignment index:  

        Consider the following alignment:  
        query residue:        1234 5678  
        query:               -ehhc-ahat-  
        target:              heh-ta-hn-g  
        alignment index:     01234567890  
                             0         1  

        For the motif 6_H_0, this target would pass. 6_H corresponds to alignment
        index 7, which is an H in the target. Note that the input POS for the motif 
        is 1-indexed.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_dali_motif_finder -h`](.github/img/aln_dali_motif_finder.png)  

# SAT aln_dali_to_pairwise_fastas
This subcommand converts the alignment fields of a DALI output file to pairwise fasta files, one per alignment. The usecase here is typically when you are trying to generate a structure-guided MSA. This means that, during the DALI  run, the 'alignments' output field must have been included.

For each pairwise alignment in the DALI output, a separate FASTA file is created with the naming pattern [query_id]xxx[target_id].fasta.aln in the specified output directory. Each FASTA file contains two sequences:  
1. The query sequence from the alignment
2. The target sequence from the alignment  

If a structure key file is provided, DALI IDs will be converted to the corresponding structure names found in the key file.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_dali_to_pairwise_fastas -h`](.github/img/aln_dali_to_pairwise_fastas.png)  

# SAT aln_filter
This subcommand filters a foldseek alignment file to keep only those alignments with a value below/above the specified value in a field (alntmscore is a common one). It also only outputs a maximum of N alignments for each query.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_filter -h`](.github/img/aln_filter.png)  

# SAT aln_merge
This subcommand is used to merge two foldseek alignment files.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_merge -h`](.github/img/aln_merge.png)  

# SAT aln_merge_clusters
This subcommand merges two cluster files by adding a higher-level clustering to a nested cluster file. 

**File1** contains the higher-level clustering (e.g., `lol_rep` -> `struc_rep`).  
**File2** contains the lower-level/nested clustering (e.g., `struc_rep` -> `seq_rep` -> `member`).  

The second column of file1 is joined to the first column of file2. File suffixes (e.g., `.pdb`, `.fasta`) are automatically removed from all values.

If a value in file2's first column is not found in file1's second column, it will be used as the higher-level cluster representative. This would be the case, for example, if you generated the higher level cluster file from an alignment and the protein didn't have any alignments.

**Example:**

File1 (higher-level clustering):
```
lol_rep    struc_rep
A          A
A          D
E          E
```

File2 (nested clustering):
```
struc_rep  seq_rep  member
A          B        B
A          B        C
D          D        D
E          E        E
E          F        F
```

Output:
```
lol_rep    struc_rep  seq_rep  member
A          A          B        B
A          A          B        C
A          D          D        D
E          E          E        E
E          E          F        F
```

Column names can be parsed from the first line of each file (default), or provided via `--file1_colnames` and `--file2_colnames` as comma-delimited lists.

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_merge_clusters -h`](.github/img/aln_merge_clusters.png)  

# SAT aln_parse_dali
This subcommand reads in a DALI alignment output file and formats it as a tab-delimited file. This script will written to the specified output file. There is also functionality to filter the alignments by zscore, alnlen, coverage, or rmsd.  
There are two main inputs:  
1) alignment_file: This is the DALI alignment file. Notably, the first output field MUST BE the 'summary' and the second output field MUST BE 'equivalences'
2) structure_key: DALI only processes files that have a 4-digit identifier. The structure key must be of format structure[delimiter]identifier, and lets you convert the identifiers back to the actual structure name. Note that the structure_key identifiers should not have the DALI segment (e.g. A, B, C...) at the end - this will be taken care of.  

The qlen field of the output is dependent on their being a self alignment in the alignment file, as then the qlen=tlen. If not present, qlen will be listed as 0.  
            
Note also the coverage is determined by alnlen/max(qlen, tlen)  

The output file is a .m8 file (e.g. tab delimited) and has the following columns: query, target, query_id, target_id, alnlen, qlen, tlen, cov, pident, rmsd, z
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_parse_dali -h`](.github/img/aln_parse_dali.png)  

# SAT aln_parse_dali_matrix
This subcommand takes in a DALI matrix file and/or a DALI dendogram files, and uses the
specified key to convert each ID to its proper name.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_parse_dali_matrix -h`](.github/img/aln_parse_dali_matrix.png)  

# SAT aln_query_uniprot
This script takes alphafold IDs (or raw uniprot IDs) and uses the Uniprot REST API to get information on the geneName and fullName (an informative name of the protein) for each ID. You can specify in which column of the infile the IDs live.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_query_uniprot -h`](.github/img/aln_query_uniprot.png)  

# SAT aln_taxa_counts
This takes in a cluster file (required columns are cluster_ID, cluster_rep, cluster_member, and cluster_count) and tallies up the taxons for each cluster. It makes a tidy file for each cluster where, for every taxon at every level, it specifies the count. The cluster file is assumed to be generated from an all-by-all alignment, perhaps with some additional merging steps. If you are also interested in adding taxonomy count information for the targets of a search of the cluster members against a separate database, you can enter an alignment file to this script. In the event an alignment file is provided, taxonIDs from the TARGET will be added to the cluster_ID of the QUERY.  
            
The output file has the following columns:  
cluster_ID, cluster_rep, cluster_count, superkingdom, level, taxon, count.  

<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_taxa_counts -h`](.github/img/aln_taxa_counts.png)  

# SAT aln_trim_target_fasta
This subcommand trims the fasta sequences of target accessions from alignment files to the longest (0) or shortest(1) target alignment length.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py aln_trim_target_fasta -h`](.github/img/aln_trim_target_fasta.png) 

# SAT plot_pae
This subcommand produces a PAE matrix plot when given a colabfold scores json file. The output file type is specified by the suffix of the out_image arguemnt.    
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py plot_pae -h`](.github/img/plot_pae.png)  

# SAT seq_chunk
Splits entries into a fasta into overlapping or non-overlapping chunks. This is helpful when you want to split up sequences that are too long to effectively use for structure prediction. This subcommand is able to generate overlapping sequences. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_chunk -h`](.github/img/seq_chunk.png)  

# SAT seq_multimerize  
 This subcommand combines input fastas to generate a multimierzed fasta containing :'s separating sequence. The cardinality of the input files can be specified to generate different kinds of homo- or hetero-complexes.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_multimerize -h`](.github/img/seq_multimerize.png)  

# SAT seq_parse_genbank   
This subcommand parses a genbank file (based on a nuclear accesion!) into an output fasta and an output table.

The output fasta will have headers with the following information:
{genome_acc}{args.delimiter}{protein_id}{args.delimiter}{locus_tag}{args.delimiter}{protein_order}
This is equivelant to the "output_name"

The output table is CSV FORMATTED, with the following columns:
{output_name},{genome_acc},{locus_tag},{protein_id},{start},{end},{strand},{protein_order},{organism_name},{protein_name}

Note that if you desire to only process a subset of genbank entires, you can provide a file with the genome accessions (no version!) that you desire.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_parse_genbank -h`](.github/img/seq_parse_genbank.png)  

# SAT seq_split_fasta
This subcommand parses a fasta file with multiple entries into multiple fasta files with one entry each.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py seq_split_fasta -h`](.github/img/seq_split_fasta.png) 

# SAT struc_clean
This subcommand takes in a structure predction (PDB file), removes all the molecules that are not amino acids (water, ions, small ligands), and outputs the cleaned protein-only PDB file.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_clean -h`](.github/img/struc_clean.png) 

# SAT struc_detect_interaction
This subcommand takes in a structure predction (PDB file) and its associated PAE file (from colabfold) that was generated with AF Multimer between two molecules. Thus, the structure prediction should have two chains, A and B. This script clusters the PAE matrix and determins if a cluster contains residues from both chains - if so, the molecules are considered to have an interaction. This script also counts the number of residues of each chain that have a C-alpha  within a specified agstrom distance from a C-alpha from the other chain.  

It is assumed the input structure has a delimiter which indicates the two members that were folded together - this delimiter can be provided.  

The output file is tab-delimited with the following columns:  
- member1  
- member2  
- interaction (True or False)  
- number of residues in chain1
- number of residues in chain2
- number of residues in chain1 that are present in cross-chain clusters  
- number of residues in chain2 that are present in corss-chain clusters  
- fraction of residues in chain 1 that are present in cross-chain clusters  
- fraction of residues in chain 2 that are present in cross-chain clusters  
- The number of residues in chain1 that have a C-alpha within distance_cutoff angstroms of a C-alpha from chain2.  
- The number of residues in chain2 that have a C-alpha within distance_cutoff angstroms of a C-alpha from chain1.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_detect_interaction -h`](.github/img/struc_detect_interaction.png) 

# SAT struc_disorder
This takes an input structure and calculates the number of residues that are considered ordered, disordered, or intermediate. A residue is considered ordered if it is in a stretch of at least n_sequential residues that have a pLDDT of >= order_cutoff. A residue is considered disordered if it is in a stretech of at least n_sequential residues <= disorder_cutoff.  

This returns an output file with the following columns:  
- basename of the input structure  
- number of ordered residues  
- number of disordered residues   
- number of intermediate residues (neither ordered or disordered)  
- total number of residues  
- there_is_a_domain: yes or no. This checks that there is at least one stretech of continuous residues that have ordered pLDDTs. The required stretch size is args.check_for_domain_len  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_disorder -h`](.github/img/struc_disorder.png) 

# SAT struc_download
This subcommand takes in a file of uniprot IDs and downloads the AF2 database pdb and pae files to the indicated directory. Furthermore, if any additional information is present in the tabular infile it will be appended to the output files - this is a good way to lable the files with information like taxonomyID, etc.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_download -h`](.github/img/struc_download.png)  

# SAT struc_extract_chains
This subcommand extracts one or more chains from an input structure, and writes them to a new pdb file. The desired chains should be input as a comma-delimited string.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_extract_chains -h`](.github/img/struc_extract_chains.png) 

# SAT struc_find_motif
Given a motif of structure [OPTIONS]xxx[OPTIONS]xx, where x indicates any amino acid and [] indicate any of the amino acids present within the brackets, this returns the match and position start/end of the motif present in the input sequence.  

The input can be a structure file, a fasta, or just a sequence. The output is tab-delimited and printed to the screen, with the columns  
- match  
- start   
- end   
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_find_motif -h`](.github/img/struc_find_motif.png) 

# SAT struc_get_contact_probability  
This script calculates the internal contact probability between two chains folded using colabfold. This can be from any model type (e.g. alphafold2_multimer_v3, or even alphafold2_ptm) as long as two chains were predicted using colabfold.  

This script parses the entire colabfold output directory. For a given model type, if different models were run (between 1-5) this script will parse every model.  

This will measure the highest contact probability in the submatrix of residues that cross chains A and B. The user can set a distance cutoff - if set to e.g. 8, this script, for every residue-residue pair between chains, will determine the probability the residue pairs are within 8 angstroms. The highest residue-residue contact probability is the result. I've annotated the functions in this script if you want to learn more about the technical details.  

For each model, there are required input files:   
1) The a3m file. This is only present once regardless of the number of models run. The critical part of this file is the first line, which lists the length of the two protein chains.  
2) The pickle file. This is the ouput pickle object from each model. This is the source of the logits from which we calculate interaction probability.  
3) The json file. This is the scores json file that lists PAE, but also at the end contains the iptm.  

The output file is tab delimited and has the following columns, with one line per model:  
- sample_name (the prefix of the input files)  
- model type (e.g. alphafold2_ptm)  
- rank (1-5, ranked based on however colabfold was ranking structures. Usually pLDDT. E.g. rank 1 is the structure with highest pLDDT)  
- model (1-5, this is the model run)  
- iptm (interaction pTM score reported by colabfold) 
- contact probability (highest cross-chain residue-residue contact probability as discussed above)   
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_get_contact_probability -h`](.github/img/struc_get_contact_probability.png)  

# SAT struc_get_domains
Extract separate domain structures from a predicted structure.  
This uses the chainsaw predicted domain boundaries to extract domains from pdb structure files. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_get_domains -h`](.github/img/struc_get_domains.png)  

# SAT struc_get_iptm  
Extract the iPTM value from the colabfold json. Appends to an output file with format {json_basename}\t{iptm}\n. 
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_get_iptm -h`](.github/img/struc_get_iptm.png)  

# SAT struc_qc
Given a structure, determines the percentage of residues that have at least the specified pLDDT. The output is returned to STDOUT!! It is tab-delimited and has the following columns:   
- structure_file (the basename of the file, including suffix)  
- number of residues  
- number of residues that pass the pLDDT threshold  
- proportion of residues that pass the pLDDT threshold (this will be a decimal between 0 and 1)
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_qc -h`](.github/img/struc_qc.png) 

# SAT struc_rebase
Simple subcommand that renumbers all residues in a structure such that the first residue is #1 and all residues are sequential (e.g. it takes out numeric gaps in residue numbers).
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_rebase -h`](.github/img/struc_rebase.png)  

# SAT struc_to_plddt
Simple subcommand that returns the average plddt of the input structure file. If --out_file is not specified, the average plddt is simply printed to the screen. If --out_file is specified, the output file will be APPENDED to with the following: [basename input structure_file]\\t[plddt]\\n
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_to_plddt -h`](.github/img/struc_to_plddt.png)  

# SAT struc_to_seq
Simple subcommand to produce the amino-acid sequence from a structure file.  

Can append the sequence to an outfile if provided, or will print to screen.
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py struc_to_seq -h`](.github/img/struc_to_seq.png)  

# SAT tab_add_taxonomy  
Adds taxonomy information to any file. There must be a column named 'taxid', that contains each taxon ID. This file will essentially add additional columns, one per desired taxonomy level, with the corresponding taxon names.  
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py tab_add_taxonomy -h`](.github/img/tab_add_taxonomy.png)  


# Planned improvements
ete3  
- Add ability to specify where the ete3 taxonomy database is downloaded.
