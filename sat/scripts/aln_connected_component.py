import random

from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def parse_alignment_file(filepath, colnames=""):
    """
    Parse an alignment file and return query-target pairs.

    If colnames is empty, the first line must start with "query" and will be
    used as the header. Otherwise, colnames should be a comma-delimited string.

    Args:
        filepath: Path to the alignment file
        colnames: Optional comma-delimited string of column names

    Returns:
        list of (query, target) tuples
    """
    pairs = []

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

    # Find query and target column indices
    if "query" not in columns:
        raise ValueError("Column names must include 'query'")
    if "target" not in columns:
        raise ValueError("Column names must include 'target'")

    query_idx = columns.index("query")
    target_idx = columns.index("target")

    # Parse data rows
    for line in lines[data_start:]:
        parts = line.split("\t")
        if len(parts) < max(query_idx, target_idx) + 1:
            continue  # Skip malformed lines

        query = parts[query_idx]
        target = parts[target_idx]
        pairs.append((query, target))

    return pairs


def find_connected_components(pairs, all_inputs=None):
    """
    Find connected components from query-target pairs.

    Uses union-find algorithm for efficient clustering.

    Args:
        pairs: list of (query, target) tuples
        all_inputs: optional set of all input members. Members not present in
                    pairs will be added as single-member clusters.

    Returns:
        list of sets, where each set is a connected component (cluster)
    """
    # Parent dictionary for union-find
    parent = {}

    def find(x):
        if x not in parent:
            parent[x] = x
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Build union-find structure
    for query, target in pairs:
        union(query, target)

    # Track all seen members
    seen_members = set(parent.keys())

    # Group by root parent
    components = {}
    for member in parent:
        root = find(member)
        if root not in components:
            components[root] = set()
        components[root].add(member)

    clusters = list(components.values())

    # Add single-member clusters for inputs not in any alignment
    if all_inputs is not None:
        not_seen = all_inputs - seen_members
        for member in not_seen:
            clusters.append({member})

    return clusters


def write_cluster_file(output_file, clusters):
    """
    Write clusters to a foldseek-style cluster file.

    Args:
        output_file: Path to output file
        clusters: list of sets, each set is a cluster
    """
    make_output_dir(output_file)

    with open(output_file, "w") as f:
        # Write header
        f.write("cluster_rep\tcluster_member\n")

        # Write data rows
        for cluster in clusters:
            # Pick a random cluster rep
            cluster_list = list(cluster)
            cluster_rep = random.choice(cluster_list)

            for member in cluster:
                f.write(f"{cluster_rep}\t{member}\n")


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_connected_component_main(args):
    """
    Generate connected component clusters from an alignment file.

    All query-target pairs that are connected (directly or transitively)
    will be placed into the same cluster.
    """

    talk_to_me("Parsing alignment file.")
    pairs = parse_alignment_file(args.alignment_file, args.colnames)
    talk_to_me(f"Found {len(pairs)} query-target pairs.")

    # Load all_inputs if provided
    all_inputs = None
    if args.all_inputs != "":
        talk_to_me("Loading all_inputs file.")
        all_inputs = set()
        with open(args.all_inputs) as f:
            for line in f:
                line = line.rstrip("\n")
                if line:
                    all_inputs.add(line)
        talk_to_me(f"Loaded {len(all_inputs)} members from all_inputs file.")

    talk_to_me("Finding connected components.")
    clusters = find_connected_components(pairs, all_inputs)

    if all_inputs is not None:
        # Count how many were added as singletons
        seen_in_pairs = set()
        for q, t in pairs:
            seen_in_pairs.add(q)
            seen_in_pairs.add(t)
        singletons_added = len(all_inputs - seen_in_pairs)
        talk_to_me(f"Added {singletons_added} members as single-member clusters.")

    talk_to_me(f"Found {len(clusters)} clusters.")

    talk_to_me(f"Writing output to {args.output_file}")
    write_cluster_file(args.output_file, clusters)

    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)

