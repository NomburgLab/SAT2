# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import pandas as pd
import os

from .utils.misc import talk_to_me, make_output_dir
from .utils.structure import pdb_to_structure_object, write_structure_subset

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

def parse_chainsaw_file(chainsaw_file):
    """
    Given text file of chainsaw results, read the results into a table and 
    return a dictionary of structure name:domain boundary key value pairs.
    If a structure has zero domains, it will not be included in the dictionary.

    - chainsaw_file: chainsaw results text file, containing chain_id, sequence_md5, nres, ndom, chopping, confidence, time_sec
    - domains_dict: dictionary of structure name: domain boundaries
    """
    chainsaw_table = pd.read_csv(chainsaw_file, sep='\t')
    
    talk_to_me("Number of structures with no domains: " + str(chainsaw_table['chopping'].isnull().sum()))
    talk_to_me("Number of structures with at least 1 domain: " + str(chainsaw_table['chopping'].notnull().sum()))

    chainsaw_table_filtered = chainsaw_table.dropna(subset=['chopping'])
    domains_dict = chainsaw_table_filtered.set_index('chain_id')['chopping'].to_dict()

    return domains_dict

def parse_domain(domain_boundaries):
    """
    Given domain_boundaries (e.g. 1-10), generate a list of consecutive numbers from the first domain residue to
    the last domain residue.

    Inputs:
    - domain_boundaries: boundaries from a domain or subdomain

    Outputs:
    - domain_residues: a list of domain residues
    """
    domain_start = int(domain_boundaries.split('-')[0])
    domain_end = int(domain_boundaries.split('-')[1])
    domain_residues = [i for i in range(domain_start, domain_end+1)]

    return domain_residues 


def format_chainsaw_domains(structure_name, domain_dict):    
    """
    Given the structure name and domain dictionary, retrieve the domain boundaries
    for the given structure. Use the domain boundaries to create a domain_residues_list. 
    Each list in domain_residues_list is some iterable/list of numbers, where each number is the position
    of one of the residues to keep.

    Input: 
    - structure_name: name of the structure
    - domain_dict: dictionary of structure name: domain boundaries
    
    Output:
    - domain_residues_list: a list of lists of domain residues corresponding to the structure_name. 
    """

    domain_boundaries = domain_dict.get(structure_name) 
    domain_list = domain_boundaries.split(',')
    domain_residues_list = []
    for domain in domain_list:
        subdomain_list = domain.split('_')
        subdomain_residues_list = []
        for subdomain in subdomain_list:
            subdomain_residues_list = subdomain_residues_list + parse_domain(subdomain)
        domain_residues_list.append(subdomain_residues_list)
           
    return domain_residues_list

def get_pdb_filename(file_path):
    """
    Give the path to a pdb file, extract and output the name of the pdb file.

    - file_path: path to pdb file
    """
    if file_path.endswith(".pdb"):
        input_file_name = os.path.splitext(os.path.basename(file_path))[0]
        return input_file_name
    else:
        talk_to_me("This is not a pdb file")
        return

def get_outfile_name(file_path, domain_start, domain_end):
    """
    Given the path to a pdb file and residue numbers for the start and end of the domain,
    generate the output file name

    Inputs:
    - file_path: path to pdb structure file
    - domain_start:  first residue of the domain
    - domain_end: last residue of the domain

    Output:
    - output_file_name: name of the extracted domain output file
    """
    input_file_name = get_pdb_filename(file_path)
    output_file_name = input_file_name + '_domain' + '_'+ str(domain_start) + '_' +  str(domain_end) + '.pdb'
    
    return output_file_name


def struc_extract_residues(pdb_file_path, domain_residues_list, min_domain_length, outfile_dir):
    """
    For the pdb file, extract and output the domains (pdb format) that meet the min_domain_length requirement.
    
    - pdb_file_path: path to pdb file
    - domain_residues_list: a list of lists of domain residues
    - min_domain_length: The length cutoff for the domains. 
                         If domain length < min_domain_length, the domain will not be written out to file.
    - outfile_dir: directory to output the pdb files
    """
    structure = pdb_to_structure_object(pdb_file_path, structure_name="structure")
    for domain_residues in domain_residues_list:
        domain_start = min(domain_residues)
        domain_end = max(domain_residues)

        if len(domain_residues) < min_domain_length:
            continue
        else:
            output_file_name = get_outfile_name(pdb_file_path, domain_start, domain_end)
            file_path = os.path.join(outfile_dir, output_file_name) 
            write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)

def struc_get_domains_main(args):
    domain_dict = parse_chainsaw_file(args.chainsaw_file_path)
    structure_name = get_pdb_filename(args.structure_file_path)
    if domain_dict.get(structure_name) is None:
        talk_to_me(f"{structure_name}PDB file is not in domain dictionary. It does not have a domain.")
    elif domain_dict.get(structure_name) != None:
        domain_residues_list = format_chainsaw_domains(structure_name=structure_name, domain_dict=domain_dict)
        make_output_dir(args.outfile_dir, is_dir=True)
        struc_extract_residues(args.structure_file_path, domain_residues_list, args.min_domain_length, args.outfile_dir)

if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
