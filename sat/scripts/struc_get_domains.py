# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import os

from .utils.misc import talk_to_me, make_output_dir, read_tsv
from .utils.structure import pdb_to_structure_object, parse_structure_inputs, write_structure_subset

# ------------------------------------------------------------------------------------ #
# Constants
# ------------------------------------------------------------------------------------ #

_STRUCTURE_EXTENSIONS = ('.pdb', '.cif', '.mmcif', '.ent')
_NULL_VALUES = frozenset({'null', 'none', 'nan', '', 'na', 'n/a', 'missing'})

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

def _strip_structure_extension(name):
    """Strip common structure file extensions so IDs match PDB basenames."""
    for ext in _STRUCTURE_EXTENSIONS:
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return name

def _build_domain_dict(rows, domain_column, id_column=None):
    """
    Build a lookup dict {structure_name: boundary_or_None} by iterating over
    rows from read_tsv (or any iterable of dicts). Each row is processed and
    discarded — only the extracted ID and boundary are kept in memory.

    Args:
        rows: iterable of dicts (e.g. read_tsv generator or a list)
        domain_column: name of the column containing domain boundaries
        id_column: name of the column containing structure IDs (default: first column)

    Returns:
        dict: {stripped_structure_name: boundary_string_or_None}
    """
    domain_dict = {}
    for row in rows:
        # Resolve id column on first row
        if id_column is None:
            id_column = next(iter(row))

        if domain_column not in row:
            raise ValueError(
                f"Domain column '{domain_column}' not found. "
                f"Available columns: {list(row.keys())}"
            )

        structure_id = _strip_structure_extension(row[id_column].strip())
        boundary = row[domain_column].strip()
        is_null = boundary.lower() in _NULL_VALUES
        domain_dict[structure_id] = None if is_null else boundary

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
                         Only one domain is accepted (e.g. 1-10 or 10-20_50-100). Underscores denote a discontinous domain.")

    domain_residues =[]
    subdomains = domain_boundary.split('_')
    for subdomain in subdomains:

        if '-' not in subdomain:
            # Single-residue segment (e.g. "20" means just residue 20)
            subdomain_start = int(subdomain)
            subdomain_end = subdomain_start
        else:
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
    - domain_boundary: a string that indicates the domain boundary (e.g. '10-20' or '10-20_50-100'),
      or 'FULL' for a single domain spanning the whole chain, or 'UNK' when no domain info is available.

    Output:
    - output_file_name: name of the extracted domain output file
    """
    input_file_name = get_pdb_filename(file_path)
    output_file_name = input_file_name + '__D' + domain_boundary +'.pdb'

    return output_file_name


def _process_one_structure(structure_file_path, domain_dict, min_domain_length, outfile_dir):
    """
    Process a single PDB structure: look up its domain boundaries and extract domains.

    Args:
        structure_file_path: path to the .pdb file
        domain_dict: {structure_name: boundary_or_None}
        min_domain_length: minimum residue count for a domain to be extracted
        outfile_dir: output directory path
    """
    structure = pdb_to_structure_object(structure_file_path, structure_name="structure")
    structure_name = get_pdb_filename(structure_file_path)
    all_pdb_residues = set(res.id[1] for res in structure.get_residues())

    if structure_name not in domain_dict:
        raise ValueError(f"{structure_name} is not found in the domain file.")

    if domain_dict[structure_name] is None:
        talk_to_me(f"{structure_name} has no domain boundary info (null). Falling back to full structure (__DUNK).")
        domain_residues = list(all_pdb_residues)
        output_file_name = get_outfile_name(structure_file_path, domain_boundary="UNK")
        file_path = os.path.join(outfile_dir, output_file_name)
        write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)
        talk_to_me(f"{structure_name}: 1 extracted domain (full chain, unknown segmentation).")
        return

    domain_boundaries = domain_dict[structure_name].split(",")
    talk_to_me(f"{structure_name} has {len(domain_boundaries)} domains.")

    ndom_extracted = 0
    for domain_boundary in domain_boundaries:

        domain_residues = parse_domain(domain_boundary)

        if len(domain_residues) < min_domain_length:
            talk_to_me(f"{structure_name} with {domain_boundary} domain boundary does not meet the minimum domain length. This domain will not be extracted.")
            continue

        ndom_extracted += 1

        # Determine output label
        is_full_chain = (set(domain_residues) == all_pdb_residues)
        if is_full_chain:
            output_file_name = get_outfile_name(structure_file_path, domain_boundary="FULL")
        else:
            output_file_name = get_outfile_name(structure_file_path, domain_boundary=domain_boundary)

        file_path = os.path.join(outfile_dir, output_file_name)
        write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)

    talk_to_me(f"{structure_name} has {ndom_extracted} extracted domains.")


def struc_get_domains_main(args):
    """
    Extract domains from PDB structures using a TSV domain segmentation file.

    - structure_file_path: path to a .pdb file or a directory of .pdb files
    - domain_file_path: path to a TSV domain segmentation file
    - colnames: comma-delimited column names (or empty to use first line as header)
    - domain_column: name of the column containing domain boundaries
    - id_column: name of the column containing structure IDs (default: first column)
    - min_domain_length: minimum domain length to extract
    - outfile_dir: directory to output the pdb files
    """
    pdb_files = parse_structure_inputs(args.structure_file_path)
    colnames = args.colnames
    id_column = args.id_column
    domain_column = args.domain_column

    make_output_dir(args.outfile_dir, is_dir=True)

    if len(pdb_files) > 1:
        talk_to_me(f"Found {len(pdb_files)} PDB files in {args.structure_file_path}")

        domain_dict = _build_domain_dict(
            read_tsv(args.domain_file_path, colnames=colnames),
            domain_column, id_column,
        )

        for pdb_file in pdb_files:
            _process_one_structure(pdb_file, domain_dict, args.min_domain_length, args.outfile_dir)
    else:
        # Single file: stream domain file, stop at first match
        structure_name = get_pdb_filename(pdb_files[0])
        found = False

        for row in read_tsv(args.domain_file_path, colnames=colnames):
            row_id_col = id_column if id_column is not None else next(iter(row))

            if domain_column not in row:
                raise ValueError(
                    f"Domain column '{domain_column}' not found. "
                    f"Available columns: {list(row.keys())}"
                )

            row_id = _strip_structure_extension(row[row_id_col].strip())
            if row_id == structure_name:
                boundary = row[domain_column].strip()
                is_null = boundary.lower() in _NULL_VALUES
                domain_dict = {structure_name: None if is_null else boundary}
                found = True
                break

        if not found:
            raise ValueError(f"{structure_name} is not found in the domain file.")

        _process_one_structure(pdb_files[0], domain_dict, args.min_domain_length, args.outfile_dir)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
