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


def format_chainsaw_domains(structure_name, domain_dict):    
    """
    Given the structure name and domain dictionary, retrieve the domain boundaries
    for the given structure. Use the domain boundaries to create a domain_residues_list. 
    Each list in domain_residues_list is some iterable/list of numbers, where each number is the position
    of one of the residues to keep.

    - structure_name: name of the structure
    - domain_dict: dictionary of structure name: domain boundaries
    - domain_residues_list: a list of lists of domain residues corresponding to the structure_name. 
    """
    
    def parse_domain(domain):
        domain_start = int(domain.split('-')[0])
        domain_end = int(domain.split('-')[1])
        domain_residues = [i for i in range(domain_start, domain_end+1)]
        return domain_residues 

    domain_boundaries = domain_dict.get(structure_name) 
    domain_list = domain_boundaries.split(',')
    domain_residues_list = []
    for domain in domain_list:
        if '_' not in domain:
            domain_residues_list.append(parse_domain(domain))
        elif '_' in domain:
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
 
def struc_extract_residues(pdb_file_path, domain_residues_list, min_domain_length, output_dir):
    """
    For the pdb file, extract and output the domains (pdb format) that meet the min_domain_length requirement.
    
    - pdb_file_path: path to pdb file
    - domain_residues_list: a list of lists of domain residues
    - min_domain_length: The length cutoff for the domains. 
                         If domain length < min_domain_length, the domain will not be written out to file.
    - output_dir: directory to output the pdb files
    """
    print("hello")
    structure = pdb_to_structure_object(pdb_file_path, structure_name="structure")
    for domain_residues in domain_residues_list:
        domain_start = min(domain_residues)
        domain_end = max(domain_residues)

        print("type for min_domain_length", min_domain_length)
        print(min_domain_length)
        if len(domain_residues) < min_domain_length:
            continue
        else:
            input_file_name = get_pdb_filename(pdb_file_path)
            output_file_name = input_file_name + '_domain' + '_'+ str(domain_start) + '_' +  str(domain_end) + '.pdb'

            file_path = os.path.join(output_dir, output_file_name) 
            write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)
    return

def struc_get_domains_main(args):
    domain_dict = parse_chainsaw_file(args.chainsaw_file_path)
    structure_name = get_pdb_filename(args.structure_file_path)
    domain_residues_list = format_chainsaw_domains(structure_name=structure_name, domain_dict=domain_dict)
    make_output_dir(args.outfile_dir, is_dir=True)
    struc_extract_residues(args.structure_file_path, domain_residues_list, args.min_domain_length, args.outfile_dir)
    return

if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
