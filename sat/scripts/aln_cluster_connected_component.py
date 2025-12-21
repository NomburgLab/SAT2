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

    # Build union-find structure - skip self-alignments as they don't connect anything
    for query, target in pairs:
        if query != target:  # Only process non-self alignments
            union(query, target)
        else:
            # Still need to add self-aligned members to parent
            find(query)

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


def analyze_alignment_connectivity(pairs):
    """
    Diagnostic function to analyze what's connecting clusters together.

    Returns information about which alignments are responsible for
    merging clusters. Useful for debugging why you might get one giant cluster.

    Args:
        pairs: list of (query, target) tuples

    Returns:
        dict with diagnostic information
    """
    # Count non-self alignments
    non_self_alignments = [(q, t) for q, t in pairs if q != t]
    self_alignments = [(q, t) for q, t in pairs if q == t]

    # Find unique members
    all_members = set()
    for q, t in pairs:
        all_members.add(q)
        all_members.add(t)

    # Build adjacency list for non-self alignments
    adjacency = {}
    for q, t in non_self_alignments:
        if q not in adjacency:
            adjacency[q] = set()
        if t not in adjacency:
            adjacency[t] = set()
        adjacency[q].add(t)
        adjacency[t].add(q)

    # Find highly connected nodes (potential "hub" proteins)
    connection_counts = {member: len(adjacency.get(member, set())) for member in all_members}
    sorted_by_connections = sorted(connection_counts.items(), key=lambda x: -x[1])

    return {
        "total_pairs": len(pairs),
        "self_alignments": len(self_alignments),
        "non_self_alignments": len(non_self_alignments),
        "unique_members": len(all_members),
        "members_with_only_self_alignment": len(all_members) - len(adjacency),
        "top_connected_members": sorted_by_connections[:10],
        "non_self_alignment_examples": non_self_alignments[:20],
    }


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
def aln_cluster_connected_component_main(args):
    """
    Generate connected component clusters from an alignment file.

    All query-target pairs that are connected (directly or transitively)
    will be placed into the same cluster.
    """

    talk_to_me("Parsing alignment file.")
    pairs = parse_alignment_file(args.alignment_file, args.colnames)
    talk_to_me(f"Found {len(pairs)} query-target pairs.")

    # Diagnostic mode
    if args.diagnose:
        talk_to_me("Running diagnostic analysis...")
        diag = analyze_alignment_connectivity(pairs)

        print("\n" + "=" * 60)
        print("DIAGNOSTIC REPORT")
        print("=" * 60)
        print(f"Total alignment pairs: {diag['total_pairs']}")
        print(f"Self-alignments (A->A): {diag['self_alignments']}")
        print(f"Non-self alignments (A->B): {diag['non_self_alignments']}")
        print(f"Unique members: {diag['unique_members']}")
        print(f"Members with ONLY self-alignment: {diag['members_with_only_self_alignment']}")
        print()
        print("Top 10 most connected members (potential hubs):")
        for member, count in diag["top_connected_members"]:
            print(f"  {member}: {count} connections")
        print()
        if diag["non_self_alignment_examples"]:
            print("First 20 non-self alignments (these create connections):")
            for q, t in diag["non_self_alignment_examples"]:
                print(f"  {q} -> {t}")
        else:
            print("No non-self alignments found!")
            print("This means all proteins should be in separate clusters.")
        print("=" * 60 + "\n")

        talk_to_me("Diagnostic complete. Continuing with clustering...")

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

    # Report largest cluster size
    if clusters:
        largest_cluster_size = max(len(c) for c in clusters)
        talk_to_me(f"Largest cluster has {largest_cluster_size} members.")

    talk_to_me(f"Writing output to {args.output_file}")
    write_cluster_file(args.output_file, clusters)

    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)

