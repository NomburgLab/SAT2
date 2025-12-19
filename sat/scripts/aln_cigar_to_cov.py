import re

from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def parse_cigar_matches(cigar_string):
    """
    Parse a CIGAR string and count the total number of matches (M operations).

    CIGAR strings consist of operations like:
    - M: match/mismatch (aligned residues)
    - I: insertion to reference
    - D: deletion from reference

    Args:
        cigar_string: A CIGAR string (e.g., "118M", "6M9I19M3I1M1I1M1I19M1D8M1D30M")

    Returns:
        int: Total count of M (match) operations
    """
    # Pattern to match numbers followed by operation codes
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")
    matches = pattern.findall(cigar_string)

    total_matches = 0
    for count_str, op in matches:
        if op == "M":
            total_matches += int(count_str)

    return total_matches


def parse_alignment_file(filepath, colnames=""):
    """
    Parse an alignment file and return column names and data rows.

    If colnames is empty, the first line must start with "query" and will be
    used as the header. Otherwise, colnames should be a comma-delimited string.

    Args:
        filepath: Path to the alignment file
        colnames: Optional comma-delimited string of column names

    Returns:
        tuple: (columns, data_rows) where columns is a list of column names
               and data_rows is a list of lists (each inner list is a row's values)
    """
    with open(filepath) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    if not lines:
        raise ValueError(f"File {filepath} is empty!")

    first_line_parts = lines[0].split("\t")

    # Determine column names and where data starts
    if colnames == "":
        # Auto-detect: first line must start with "query"
        if first_line_parts[0] == "query":
            columns = first_line_parts
            data_start = 1
            talk_to_me(f"Auto-detected column names from header: {columns}")
        else:
            msg = (
                "colnames not provided and first line does not start with 'query'. "
                "Either provide column names via --colnames or ensure the alignment "
                "file has a header row starting with 'query'."
            )
            raise ValueError(msg)
    else:
        # Use provided colnames
        columns = [c.strip() for c in colnames.split(",")]
        talk_to_me(f"Using provided column names: {columns}")

        # Check if first line is a header that matches our colnames
        if first_line_parts[0] == "query":
            talk_to_me("First line appears to be a header, skipping it.")
            data_start = 1
        else:
            data_start = 0

    # Validate required columns exist
    required_cols = ["qlen", "tlen", "cigar"]
    for col in required_cols:
        if col not in columns:
            raise ValueError(
                f"Required column '{col}' not found in alignment fields: {columns}"
            )

    # Parse data rows
    data_rows = []
    for line in lines[data_start:]:
        parts = line.split("\t")
        data_rows.append(parts)

    return columns, data_rows


def add_cigar_coverage_columns(columns, data_rows):
    """
    Add cigar_qcov and cigar_tcov columns to the data.

    The new columns are inserted directly after the 'cigar' column.

    Args:
        columns: List of column names
        data_rows: List of data rows (each row is a list of values)

    Returns:
        tuple: (new_columns, new_data_rows) with cigar_qcov and cigar_tcov added
    """
    qlen_idx = columns.index("qlen")
    tlen_idx = columns.index("tlen")
    cigar_idx = columns.index("cigar")

    # Insert new column names after cigar
    new_columns = (
        columns[: cigar_idx + 1]
        + ["cigar_qcov", "cigar_tcov"]
        + columns[cigar_idx + 1 :]
    )

    new_data_rows = []

    for row in data_rows:
        try:
            qlen = int(row[qlen_idx])
            tlen = int(row[tlen_idx])
            cigar = row[cigar_idx]
        except (ValueError, IndexError) as e:
            talk_to_me(f"Warning: Skipping row due to parsing error: {e}")
            continue

        # Calculate coverage (handle zero length)
        num_matches = parse_cigar_matches(cigar)

        if qlen > 0:
            qcov = num_matches / qlen
        else:
            qcov = 0.0

        if tlen > 0:
            tcov = num_matches / tlen
        else:
            tcov = 0.0

        # Build new row with coverage columns inserted after cigar
        new_row = (
            row[: cigar_idx + 1]
            + [f"{qcov:.3f}", f"{tcov:.3f}"]
            + row[cigar_idx + 1 :]
        )
        new_data_rows.append(new_row)

    return new_columns, new_data_rows


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_cigar_to_cov_main(args):
    """
    Add CIGAR-derived query and target coverage columns to an alignment file.

    This script:
    1. Parses the alignment file
    2. Extracts CIGAR strings and counts M (match) operations
    3. Calculates cigar_qcov = #Ms/qlen and cigar_tcov = #Ms/tlen
    4. Inserts cigar_qcov and cigar_tcov columns after the cigar column
    """
    talk_to_me("Parsing alignment file.")
    columns, data_rows = parse_alignment_file(
        args.alignment_file, args.colnames
    )
    talk_to_me(f"Found {len(data_rows)} alignment rows.")

    talk_to_me("Adding CIGAR coverage columns.")
    new_columns, new_data_rows = add_cigar_coverage_columns(columns, data_rows)
    talk_to_me(f"Processed {len(new_data_rows)} rows.")

    talk_to_me(f"Writing output to {args.output_file}")
    make_output_dir(args.output_file)

    with open(args.output_file, "w") as f:
        # Write header
        f.write("\t".join(new_columns) + "\n")
        # Write data
        for row in new_data_rows:
            f.write("\t".join(row) + "\n")

    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
