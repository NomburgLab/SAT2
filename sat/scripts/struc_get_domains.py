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
    Given text file of chainsaw results, read the results into a table and return a nested dictionary.
    The dictionary key is the structure name (i.e. chain_id in the chainsaw file). 
    Each structure name maps to an inner dictionary containing key-value pairs that
    represent the characteristics of the structure file from the chainsaw results.

    Inputs:
    - chainsaw_file: chainsaw results text file, containing chain_id, sequence_md5, nres, ndom, chopping, confidence, time_sec

    Outputs:
    - domain_dict: a nested dictionary of each structure's chainsaw result. 
                    {stucture_name: {chain_id:value, sequence_md5:value, nres:value, ndom:value, chopping:value, confidence:value, time_sec:value}}
                    Domain boundaries are found by using the "chopping" key. 
                    Structures with no domains have their "chopping" key filled with '1-nres' (from one to the total number of residues).
    """
    chainsaw_table = pd.read_csv(chainsaw_file, sep='\t')

    chainsaw_table['chopping'] = chainsaw_table['chopping'].fillna(chainsaw_table['nres'].apply(lambda nres: f'1-{nres}'))

    domain_dict = chainsaw_table.set_index('chain_id').to_dict(orient='index')

    return domain_dict

def parse_domain(domain_boundary):
    """
    Given a domain_boundary (e.g. 1-10 or 10-20_50-100), generate a list of consecutive numbers from the first domain residue to
    the last domain residue. Discontinuous domains are separated by underscores.

    Inputs:
    - domain_boundaries: boundaries from a domain

    Outputs:
    - domain_residues: a list of domain residues
    """
    if "," in domain_boundary:
        raise ValueError(f"{domain_boundary} domain boundary has a comma, which indicates more than one domain. \
                         Only one domain is accepted (e.g. 1-10 or 10-20_50-100).")
    
    domain_residues =[]
    subdomains = domain_boundary.split('_')
    for subdomain in subdomains:
        subdomain_start = int(subdomain.split('-')[0])
        subdomain_end = int(subdomain.split('-')[1])
        subdomain_residues = [i for i in range(subdomain_start, subdomain_end+1)]
        domain_residues.extend(subdomain_residues)

    return domain_residues 


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

def get_outfile_name(file_path, domain_boundary):
    """
    Given the path to a pdb file and domain boundary of a structure,
    generate the output file name

    Inputs:
    - file_path: path to pdb structure file
    - domain_boundary: a string that indicates the domain boundary (e.g. '10-20' or '10-20_50-100')
      Discontinuous domains are separated by underscores.

    Output:
    - output_file_name: name of the extracted domain output file
    """
    input_file_name = get_pdb_filename(file_path)
    output_file_name = input_file_name + '__D' + domain_boundary +'.pdb'
    
    return output_file_name

def struc_get_domains_main(args):
    """
    For the pdb file, extract and output the domains (pdb format) that meet the min_domain_length requirement.
    
    - structure_file_path: path to pdb file
    - chainsaw_file_path: path to chainsaw text file
    - min_domain_length: The length cutoff for the domains. 
                         If domain length < min_domain_length, the domain will not be written out to file.
    - outfile_dir: directory to output the pdb files
    """
    structure = pdb_to_structure_object(args.structure_file_path, structure_name="structure")
    structure_name = get_pdb_filename(args.structure_file_path)
    domain_dict = parse_chainsaw_file(args.chainsaw_file_path)

    #for a structure, create a list of domains from 'chopping' key
    domain_boundaries = domain_dict[structure_name]['chopping'].split(",")
    talk_to_me(f"{structure_name} has {len(domain_boundaries)} domains.")

    #check if structure exists
    if structure_name not in domain_dict:
        raise ValueError(f"{structure_name} is not found in the chainsaw file.")
    
    make_output_dir(args.outfile_dir, is_dir=True)
    ndom_extracted = 0

    #iterate through each domain of a structure and extract it only if it meets min length requirement
    for domain_boundary in domain_boundaries:
        domain_residues = parse_domain(domain_boundary)
        if len(domain_residues) < args.min_domain_length:
            talk_to_me(f"{structure_name} with {domain_boundary} domain boundary does not meet the minimum domain length. This domain will not be extracted.")
            continue
        else:
            ndom_extracted +=1

            #determine the domain boundary for the output file name
            if domain_boundary == "1-"+ str(domain_dict[structure_name]['nres']):
                output_file_name = get_outfile_name(args.pdb_file_path, domain_boundary="FULL")
            else:
                output_file_name = get_outfile_name(args.pdb_file_path, domain_boundary=domain_boundary)
            
            file_path = os.path.join(args.outfile_dir, output_file_name) 

            #write the extracted domain to file
            write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)
    

if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
