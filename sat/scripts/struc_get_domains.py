# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import re
import json
import pandas as pd
import os

from .utils.misc import talk_to_me, make_output_dir
from .utils.structure import pdb_to_structure_object, write_structure_subset

# ------------------------------------------------------------------------------------ #
# Constants
# ------------------------------------------------------------------------------------ #

DOMAIN_PATTERN = re.compile(r'^\d+-\d+([_,]\d+-\d+)*$')
_NULL_LIKE_COLUMN_NAMES = frozenset({'null', 'none', 'nan', '', 'na', 'n/a', 'missing'})

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

def detect_domain_column(df):
    """
    Given a DataFrame, return the name of the first column whose non-null values
    predominantly match the domain boundary pattern (e.g. '1-50', '3-10,15-45',
    '25-30_35-40,110-120'). Raises ValueError if no such column is found.
    """
    for col in df.columns:
        non_null = df[col].dropna().astype(str).str.strip()
        if non_null.empty:
            continue
        match_rate = non_null.apply(lambda v: bool(DOMAIN_PATTERN.match(v))).mean()
        if match_rate > 0.5:
            return col
    raise ValueError(
        "No domain boundary column detected in the file. "
        "Expected a column whose values match patterns like '1-50', '3-10,15-45', or '25-30_35-40'."
    )


def parse_domain_file(file_path, id_column=None):
    """
    Generic parser for any file containing domain segmentation data.
    Supports TSV, CSV (with or without header), and JSON (flat dict or list of records).

    The file may contain entries for a single structure or a concatenation of many
    structures (e.g. multiple tool outputs merged). All rows are loaded; the caller
    looks up the relevant structure by name.

    Inputs:
    - file_path: path to the domain segmentation file (.json, .tsv, .csv, or any
                 tab/comma-delimited text file)
    - id_column: name or 0-based integer index of the column containing structure IDs.
                 Defaults to the first column (index 0) if not provided.

    Outputs:
    - domain_dict: {structure_name: chopping_string_or_None}
                   None means the domain boundary was absent or null for that structure.
    """
    if file_path.endswith('.json'):
        return _parse_json_domain_file(file_path, id_column)
    else:
        return _parse_tabular_domain_file(file_path, id_column)


def _column_names_look_like_data(df):
    """
    Return True if any column name resembles data rather than a real header label.
    Detects when a headerless file was accidentally read with header=0: the 'header'
    row (which is actually data) would have null-like or domain-pattern column names.
    """
    for col in df.columns:
        col_str = str(col).strip().lower()
        if col_str in _NULL_LIKE_COLUMN_NAMES:
            return True
        if bool(DOMAIN_PATTERN.match(str(col).strip())):
            return True
    return False


def _parse_tabular_domain_file(file_path, id_column=None):
    """
    Parse a TSV or CSV domain file, auto-detecting the separator and whether
    a header row is present.
    """
    df = None
    for sep in ['\t', ',']:
        for header in [0, None]:
            try:
                candidate = pd.read_csv(file_path, sep=sep, header=header, dtype=str)
                if header == 0 and _column_names_look_like_data(candidate):
                    continue  # first row looks like data, not a real header
                detect_domain_column(candidate)
                df = candidate
                break
            except (ValueError, pd.errors.ParserError):
                continue
        if df is not None:
            break

    if df is None:
        raise ValueError(
            f"Could not detect a domain boundary column in '{file_path}'. "
            "Check that the file contains a column with values like '1-50' or '3-10,15-45'."
        )

    return _extract_domain_dict_from_df(df, id_column)


def _parse_json_domain_file(file_path, id_column=None):
    """
    Parse a JSON domain file. Supports:
    - Flat dict: {"structure_name": "boundary_string", ...}
    - List of records: [{"id": "...", "domains": "..."}, ...]
    """
    with open(file_path) as f:
        data = json.load(f)

    if isinstance(data, dict):
        non_null_values = [v for v in data.values() if v is not None]
        # Accept flat dict if all non-null values match the domain pattern (or all are null)
        if not non_null_values or all(
            isinstance(v, str) and bool(DOMAIN_PATTERN.match(v.strip()))
            for v in non_null_values
        ):
            return {_strip_structure_extension(k): (v if v is not None else None) for k, v in data.items()}
        raise ValueError(
            f"JSON file '{file_path}' is a dict but values don't match the domain boundary pattern. "
            "Expected either a flat {{name: boundary}} dict or a list of records."
        )

    if isinstance(data, list):
        df = pd.DataFrame(data).astype(str)
        return _extract_domain_dict_from_df(df, id_column)

    raise ValueError(
        f"JSON file '{file_path}' must be either a dict or a list of records."
    )


_STRUCTURE_EXTENSIONS = ('.pdb', '.cif', '.mmcif', '.ent')


def _strip_structure_extension(name):
    """Strip common structure file extensions so IDs match PDB basenames."""
    for ext in _STRUCTURE_EXTENSIONS:
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return name


def _extract_domain_dict_from_df(df, id_column=None):
    """
    Extract {structure_name: boundary_or_None} from a DataFrame using auto-detected columns.
    Structure IDs are normalized by stripping common file extensions (.pdb, .cif, etc.)
    so that tools like Merizo (which use filenames as IDs) are handled correctly.
    """
    domain_col = detect_domain_column(df)
    id_col = _resolve_id_column(df, id_column)

    domain_dict = {}
    for _, row in df.iterrows():
        structure_id = _strip_structure_extension(str(row[id_col]).strip())
        boundary = row[domain_col]
        is_null = pd.isna(boundary) or str(boundary).strip().lower() in ('nan', 'none', 'null', '')
        domain_dict[structure_id] = None if is_null else str(boundary).strip()

    return domain_dict


def _resolve_id_column(df, id_column):
    """
    Given a DataFrame and an optional id_column hint (name or int index),
    return the actual column label to use. Defaults to the first column.
    """
    if id_column is None:
        return df.columns[0]
    if isinstance(id_column, int):
        return df.columns[id_column]
    if id_column in df.columns:
        return id_column
    # Try interpreting as an integer string
    try:
        return df.columns[int(id_column)]
    except (ValueError, IndexError):
        raise ValueError(
            f"id_column '{id_column}' not found in the file columns: {list(df.columns)}"
        )


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
            raise ValueError(f"{domain_boundary} domain boundary must have a hyphen to separate subdomains (e.g. 1-10).")

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

def struc_get_domains_main(args):
    """
    For the pdb file, extract and output the domains (pdb format) that meet the min_domain_length requirement.

    - structure_file_path: path to pdb file
    - domain_file_path: path to a domain segmentation file (TSV, CSV, or JSON) containing
                        domain boundaries in the convention '1-50', '3-10,15-45', '25-30_35-40', etc.
                        Can be a single structure's output or a concatenation of many.
    - id_column: name or 0-based index of the column containing structure IDs (default: first column)
    - min_domain_length: The length cutoff for the domains.
                         If domain length < min_domain_length, the domain will not be written out to file.
    - outfile_dir: directory to output the pdb files
    """
    structure = pdb_to_structure_object(args.structure_file_path, structure_name="structure")
    structure_name = get_pdb_filename(args.structure_file_path)
    all_pdb_residues = set(res.id[1] for res in structure.get_residues())

    id_column = getattr(args, 'id_column', None)
    domain_dict = parse_domain_file(args.domain_file_path, id_column=id_column)

    make_output_dir(args.outfile_dir, is_dir=True)
    ndom_extracted = 0

    # Determine domain boundaries for this structure
    if structure_name not in domain_dict or domain_dict[structure_name] is None:
        if structure_name not in domain_dict:
            talk_to_me(f"{structure_name} not found in the domain file. Falling back to full structure (__DUNK).")
        else:
            talk_to_me(f"{structure_name} has no domain boundary info (null). Falling back to full structure (__DUNK).")

        domain_residues = list(all_pdb_residues)
        output_file_name = get_outfile_name(args.structure_file_path, domain_boundary="UNK")
        file_path = os.path.join(args.outfile_dir, output_file_name)
        write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)
        talk_to_me(f"{structure_name}: 1 extracted domain (full chain, unknown segmentation).")
        return

    domain_boundaries = domain_dict[structure_name].split(",")
    talk_to_me(f"{structure_name} has {len(domain_boundaries)} domains.")

    for domain_boundary in domain_boundaries:

        domain_residues = parse_domain(domain_boundary)

        if len(domain_residues) < args.min_domain_length:
            talk_to_me(f"{structure_name} with {domain_boundary} domain boundary does not meet the minimum domain length. This domain will not be extracted.")
            continue

        ndom_extracted += 1

        # Determine output label
        is_full_chain = (set(domain_residues) == all_pdb_residues)
        if is_full_chain:
            output_file_name = get_outfile_name(args.structure_file_path, domain_boundary="FULL")
        else:
            output_file_name = get_outfile_name(args.structure_file_path, domain_boundary=domain_boundary)

        file_path = os.path.join(args.outfile_dir, output_file_name)
        write_structure_subset(structure, residues_to_keep=domain_residues, outfile=file_path)

    talk_to_me(f"{structure_name} has {ndom_extracted} extracted domains.")

if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
