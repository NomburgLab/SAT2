"""
Greedy set cover clustering, similar to foldseek/mmseqs cluster mode 0.

This clustering method:
1. Builds an adjacency list from alignment pairs
2. Sorts members by number of alignments (descending)
3. Picks the top member as cluster rep, assigns all its alignments to its cluster
4. Removes assigned members from consideration
5. Repeats until all members are clustered
6. Optionally performs a reassignment step where members are moved to a different
   cluster if they have a better alignment score to that cluster's representative.
"""

from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def parse_alignment_file(filepath, colnames="", score_column="alntmscore"):
    """
    Parse an alignment file and return alignment data with scores.

    Args:
        filepath: Path to the alignment file
        colnames: Optional comma-delimited string of column names
        score_column: Name of the column containing alignment scores

    Returns:
        tuple: (list of (query, target, score) tuples, set of all members)
    """
    alignments = []
    all_members = set()

    with open(filepath) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    if not lines:
        raise ValueError(f"File {filepath} is empty!")

    first_line_parts = lines[0].split("\t")

    # Determine column names and where data starts
    if colnames == "":
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
        columns = [c.strip() for c in colnames.split(",")]
        talk_to_me(f"Using provided column names: {columns}")
        if first_line_parts[0] == "query":
            talk_to_me("First line appears to be a header, skipping it.")
            data_start = 1
        else:
            data_start = 0

    # Find column indices
    if "query" not in columns:
        raise ValueError("Column names must include 'query'")
    if "target" not in columns:
        raise ValueError("Column names must include 'target'")
    if score_column not in columns:
        raise ValueError(f"Column names must include '{score_column}'")

    query_idx = columns.index("query")
    target_idx = columns.index("target")
    score_idx = columns.index(score_column)

    # Parse data rows
    for line in lines[data_start:]:
        parts = line.split("\t")
        if len(parts) < max(query_idx, target_idx, score_idx) + 1:
            continue

        query = parts[query_idx]
        target = parts[target_idx]
        try:
            score = float(parts[score_idx])
        except ValueError:
            continue

        all_members.add(query)
        all_members.add(target)

        # Skip self-alignments for clustering purposes
        if query != target:
            alignments.append((query, target, score))

    return alignments, all_members


def build_adjacency_with_scores(alignments):
    """
    Build adjacency dict with best scores for each pair.

    Args:
        alignments: list of (query, target, score) tuples

    Returns:
        dict: {member: {neighbor: best_score, ...}, ...}
    """
    adjacency = {}

    for query, target, score in alignments:
        # Add query -> target
        if query not in adjacency:
            adjacency[query] = {}
        if target not in adjacency[query] or score > adjacency[query][target]:
            adjacency[query][target] = score

        # Add target -> query (symmetric)
        if target not in adjacency:
            adjacency[target] = {}
        if query not in adjacency[target] or score > adjacency[target][query]:
            adjacency[target][query] = score

    return adjacency


def greedy_set_cover_clustering(adjacency, all_members):
    """
    Perform greedy set cover clustering.

    Args:
        adjacency: dict mapping member -> {neighbor: score, ...}
        all_members: set of all members to cluster

    Returns:
        dict: {cluster_rep: set of members, ...}
    """
    clusters = {}
    remaining = all_members.copy()

    while remaining:
        # Find member with most connections to remaining members
        best_rep = None
        best_coverage = set()

        for member in remaining:
            # Coverage = neighbors that are still remaining + self
            neighbors = set(adjacency.get(member, {}).keys())
            coverage = (neighbors & remaining) | {member}

            if len(coverage) > len(best_coverage):
                best_rep = member
                best_coverage = coverage

        # Create cluster
        clusters[best_rep] = best_coverage
        remaining -= best_coverage

    return clusters


def reassign_to_better_reps(clusters, adjacency):
    """
    Reassign cluster members to better-scoring representatives.

    For each non-rep member, check if there's another cluster rep
    with a higher alignment score. If so, move the member.

    Args:
        clusters: dict {rep: set of members}
        adjacency: dict {member: {neighbor: score}}

    Returns:
        dict: updated clusters after reassignment
    """
    # Get all reps
    reps = set(clusters.keys())

    # Track reassignments
    reassignments = 0

    # Build member -> current rep mapping
    member_to_rep = {}
    for rep, members in clusters.items():
        for member in members:
            member_to_rep[member] = rep

    # Check each non-rep member
    for member, current_rep in list(member_to_rep.items()):
        if member in reps:
            continue  # Don't reassign reps

        # Get score to current rep
        current_score = adjacency.get(member, {}).get(current_rep, 0)

        # Check all other reps
        best_rep = current_rep
        best_score = current_score

        for other_rep in reps:
            if other_rep == current_rep:
                continue
            score = adjacency.get(member, {}).get(other_rep, 0)
            if score > best_score:
                best_rep = other_rep
                best_score = score

        # Reassign if better rep found
        if best_rep != current_rep:
            clusters[current_rep].remove(member)
            clusters[best_rep].add(member)
            member_to_rep[member] = best_rep
            reassignments += 1

    return clusters, reassignments


def write_cluster_file(output_file, clusters):
    """
    Write clusters to a foldseek-style cluster file.

    Args:
        output_file: Path to output file
        clusters: dict {rep: set of members}
    """
    make_output_dir(output_file)

    with open(output_file, "w") as f:
        f.write("cluster_rep\tcluster_member\n")

        for rep, members in clusters.items():
            for member in members:
                f.write(f"{rep}\t{member}\n")


# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def aln_greedy_cluster_main(args):
    """
    Perform greedy set cover clustering on an alignment file.
    """

    talk_to_me("Parsing alignment file.")
    alignments, all_members = parse_alignment_file(
        args.alignment_file, args.colnames, args.score_column
    )
    talk_to_me(f"Found {len(alignments)} non-self alignments.")
    talk_to_me(f"Found {len(all_members)} unique members.")

    # Add all_inputs if provided
    if args.all_inputs != "":
        talk_to_me("Loading all_inputs file.")
        with open(args.all_inputs) as f:
            for line in f:
                line = line.rstrip("\n")
                if line:
                    all_members.add(line)
        talk_to_me(f"Total members after adding all_inputs: {len(all_members)}")

    talk_to_me("Building adjacency structure with scores.")
    adjacency = build_adjacency_with_scores(alignments)

    talk_to_me("Performing greedy set cover clustering.")
    clusters = greedy_set_cover_clustering(adjacency, all_members)
    talk_to_me(f"Initial clustering: {len(clusters)} clusters.")

    # Report largest cluster
    largest = max(len(m) for m in clusters.values())
    talk_to_me(f"Largest cluster has {largest} members.")

    # Reassignment step
    if not args.no_reassign:
        talk_to_me("Performing reassignment step.")
        clusters, reassignments = reassign_to_better_reps(clusters, adjacency)
        talk_to_me(f"Reassigned {reassignments} members to better representatives.")

        # Report largest cluster after reassignment
        largest = max(len(m) for m in clusters.values())
        talk_to_me(f"After reassignment, largest cluster has {largest} members.")

    talk_to_me(f"Writing output to {args.output_file}")
    write_cluster_file(args.output_file, clusters)

    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)

