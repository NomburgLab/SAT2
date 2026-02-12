# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import os

from .utils.misc import talk_to_me, make_output_dir

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #


def parse_alignment_file(filepath, colnames=""):
    """
    Parse an alignment file and return column names, column indices, and data rows.

    If colnames is empty, the first line must start with "query" and will be
    used as the header. Otherwise, colnames should be a comma-delimited string.

    Args:
        filepath: Path to the alignment file
        colnames: Optional comma-delimited string of column names

    Returns:
        list of dicts with keys: query, target, qaln, taln
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
    required_cols = ["query", "target", "qaln", "taln"]
    for col in required_cols:
        if col not in columns:
            raise ValueError(
                f"Required column '{col}' not found in column names: {columns}"
            )

    query_idx = columns.index("query")
    target_idx = columns.index("target")
    qaln_idx = columns.index("qaln")
    taln_idx = columns.index("taln")

    # Parse data rows
    records = []
    for line in lines[data_start:]:
        parts = line.split("\t")
        if len(parts) < max(query_idx, target_idx, qaln_idx, taln_idx) + 1:
            continue
        records.append(
            {
                "query": parts[query_idx],
                "target": parts[target_idx],
                "qaln": parts[qaln_idx],
                "taln": parts[taln_idx],
            }
        )

    return records


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_create_pairwise_aln_fasta_main(args):
    talk_to_me("Reading alignment file")
    records = parse_alignment_file(args.alignment_file, args.colnames)
    talk_to_me(f"Found {len(records)} alignments")

    # Track stats
    written = 0
    skipped_self = 0

    for record in records:
        query = record["query"]
        target = record["target"]
        qaln = record["qaln"]
        taln = record["taln"]

        # Skip self-alignments unless --include_self is set (default: exclude)
        if not args.include_self and query == target:
            skipped_self += 1
            continue

        # Create query-specific output directory
        query_dir = os.path.join(args.output_dir, query)
        make_output_dir(query_dir, is_dir=True)

        # Write pairwise FASTA file
        outfile = os.path.join(query_dir, f"{query}xxx{target}.fasta.aln")
        with open(outfile, "w") as f:
            f.write(f">{query}\n")
            f.write(f"{qaln}\n")
            f.write(f">{target}\n")
            f.write(f"{taln}\n")

        written += 1

    talk_to_me(f"Wrote {written} pairwise FASTA files")
    if skipped_self > 0:
        talk_to_me(f"Skipped {skipped_self} self-alignments")
    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
