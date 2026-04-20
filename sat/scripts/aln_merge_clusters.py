import os
import re

from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def remove_file_suffix(value):
    """
    Remove common file extensions from a value.
    Handles extensions like .pdb, .fasta, .fa, .cif, .mmcif, .ent, .gz, etc.
    Also handles double extensions like .pdb.gz
    """
    # Common bioinformatics file extensions
    extensions = [
        r"\.pdb\.gz$",
        r"\.cif\.gz$",
        r"\.fasta\.gz$",
        r"\.fa\.gz$",
        r"\.pdb$",
        r"\.cif$",
        r"\.mmcif$",
        r"\.ent$",
        r"\.fasta$",
        r"\.fa$",
        r"\.faa$",
        r"\.fna$",
        r"\.ffn$",
        r"\.frn$",
    ]

    result = value
    for ext in extensions:
        result = re.sub(ext, "", result, flags=re.IGNORECASE)
    return result


def parse_cluster_file(filepath, colnames=None):
    """
    Parse a cluster file (simple 2-column or nested multi-column).

    Args:
        filepath: Path to the cluster file
        colnames: Optional list of column names. If None, assumes first line is header.

    Returns:
        tuple: (colnames list, list of row dicts)
    """
    rows = []

    with open(filepath) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    if not lines:
        raise ValueError(f"File {filepath} is empty!")

    # Parse column names
    first_line_parts = lines[0].split("\t")

    if colnames is None or colnames == "":
        # Assume first line is header
        colnames = first_line_parts
        talk_to_me(f"Using header from file: {colnames}")
        data_start = 1
    else:
        # Use provided colnames
        if isinstance(colnames, str):
            colnames = [c.strip() for c in colnames.split(",")]
        data_start = 0

        # Check if first line looks like a header (matches our colnames)
        if first_line_parts == colnames:
            talk_to_me("First line matches provided column names, skipping as header.")
            data_start = 1

    # Parse data rows
    for i, line in enumerate(lines[data_start:], start=data_start + 1):
        parts = line.split("\t")
        if len(parts) != len(colnames):
            raise ValueError(
                f"Line {i} has {len(parts)} columns but expected {len(colnames)}. "
                f"Line: {line}"
            )

        # Remove file suffixes from all values
        row = {col: remove_file_suffix(val) for col, val in zip(colnames, parts)}
        rows.append(row)

    return colnames, rows


def merge_cluster_files(file1_colnames, file1_rows, file2_colnames, file2_rows):
    """
    Merge two cluster files by joining file1's second column to file2's first column.

    Args:
        file1_colnames: Column names from file1 (higher-level clustering)
        file1_rows: Row dicts from file1
        file2_colnames: Column names from file2 (lower-level/nested clustering)
        file2_rows: Row dicts from file2

    Returns:
        tuple: (output_colnames, output_rows, unmatched_count)
            - unmatched_count: number of file2 rows where the first column value
              was not found in file1's second column (these get "X" as higher-level value)
    """
    # Get the join columns
    file1_join_col = file1_colnames[1]  # Second column of file1
    file2_join_col = file2_colnames[0]  # First column of file2

    # Build lookup from file1: maps file1's second column value to file1's first column value
    # Also track all other columns from file1 if there are more than 2
    file1_lookup = {}
    file1_extra_cols = file1_colnames[2:] if len(file1_colnames) > 2 else []

    for row in file1_rows:
        key = row[file1_join_col]
        if key in file1_lookup:
            # Check consistency - same key should map to same higher-level value
            existing = file1_lookup[key]
            if existing[file1_colnames[0]] != row[file1_colnames[0]]:
                raise ValueError(
                    f"Inconsistent mapping in file1: '{key}' maps to both "
                    f"'{existing[file1_colnames[0]]}' and '{row[file1_colnames[0]]}'"
                )
        else:
            file1_lookup[key] = row

    # Build output column names:
    # file1's first column + any extra file1 columns + all file2 columns
    output_colnames = [file1_colnames[0]] + file1_extra_cols + file2_colnames

    # Build output rows
    output_rows = []
    unmatched_file2_count = 0
    matched_file1_keys = set()

    for row2 in file2_rows:
        join_value = row2[file2_join_col]

        # Build output row
        out_row = {}

        if join_value in file1_lookup:
            row1 = file1_lookup[join_value]
            matched_file1_keys.add(join_value)
            # Add file1's first column (the new highest level)
            out_row[file1_colnames[0]] = row1[file1_colnames[0]]
            # Add any extra columns from file1
            for col in file1_extra_cols:
                out_row[col] = row1[col]
        else:
            # Value not found in file1 - use file2's first column value as its own
            # higher-level representative (it becomes its own cluster)
            unmatched_file2_count += 1
            out_row[file1_colnames[0]] = join_value
            for col in file1_extra_cols:
                out_row[col] = join_value

        # Add all columns from file2
        for col in file2_colnames:
            out_row[col] = row2[col]

        output_rows.append(out_row)

    # Add file1 entries that had no match in file2
    unmatched_file1_count = 0
    for key, row1 in file1_lookup.items():
        if key not in matched_file1_keys:
            unmatched_file1_count += 1
            out_row = {}
            out_row[file1_colnames[0]] = row1[file1_colnames[0]]
            for col in file1_extra_cols:
                out_row[col] = row1[col]
            # The join column gets the file1 value; remaining file2 columns get "X"
            out_row[file2_join_col] = key
            for col in file2_colnames[1:]:
                out_row[col] = "X"
            output_rows.append(out_row)

    return output_colnames, output_rows, unmatched_file2_count, unmatched_file1_count


def write_output(output_file, colnames, rows):
    """
    Write the merged cluster file.
    """
    make_output_dir(output_file)

    with open(output_file, "w") as f:
        # Write header
        f.write("\t".join(colnames) + "\n")
        # Write data rows
        for row in rows:
            values = [row[col] for col in colnames]
            f.write("\t".join(values) + "\n")


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_merge_clusters_main(args):
    """
    Merge two cluster files by adding a higher-level clustering to a nested cluster file.

    File1 contains the higher-level clustering (e.g., lol_rep -> struc_rep).
    File2 contains the lower-level/nested clustering (e.g., struc_rep -> seq_rep -> member).

    The second column of file1 is joined to the first column of file2.
    """

    talk_to_me("Parsing file1 (higher-level clustering).")
    file1_colnames, file1_rows = parse_cluster_file(args.file1, args.file1_colnames)

    talk_to_me("Parsing file2 (lower-level/nested clustering).")
    file2_colnames, file2_rows = parse_cluster_file(args.file2, args.file2_colnames)

    talk_to_me("Merging cluster files.")
    output_colnames, output_rows, unmatched_file2_count, unmatched_file1_count = (
        merge_cluster_files(
            file1_colnames, file1_rows, file2_colnames, file2_rows
        )
    )

    if unmatched_file2_count > 0:
        talk_to_me(
            f"NOTE: {unmatched_file2_count} rows in file2 had values in the first column "
            f"that were not found in file1's second column. These rows use their own "
            f"first column value as the higher-level cluster representative."
        )

    if unmatched_file1_count > 0:
        talk_to_me(
            f"NOTE: {unmatched_file1_count} entries in file1's second column had no "
            f"match in file2's first column. These were added to the output with 'X' "
            f"for the missing file2 columns."
        )

    # Sort output by the first column (highest-level rep) to keep clusters together
    talk_to_me("Sorting output by highest-level cluster representative.")
    first_col = output_colnames[0]
    output_rows.sort(key=lambda row: row[first_col])

    talk_to_me(f"Writing output with columns: {output_colnames}")
    write_output(args.output_file, output_colnames, output_rows)

    talk_to_me(f"Done! Wrote {len(output_rows)} rows to {args.output_file}")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
