import random
from collections import defaultdict

from .utils.misc import talk_to_me, make_output_dir, read_tsv


def aln_sample_main(args):
    talk_to_me("Reading alignments.")

    # Group rows by query column
    groups = defaultdict(list)
    columns = None

    rows_iter = read_tsv(args.alignment_file, args.alignment_fields)
    first_row = next(rows_iter, None)

    if first_row is None:
        # File had only a header (or was empty — but read_tsv raises on empty)
        talk_to_me("No data rows found. Writing header only.")
        make_output_dir(args.output_file)
        with open(args.alignment_file) as infile:
            header_line = infile.readline().rstrip("\n")
        with open(args.output_file, "w") as outfile:
            outfile.write(header_line + "\n")
        return

    columns = list(first_row.keys())

    if args.query_column not in columns:
        raise ValueError(
            f"Query column '{args.query_column}' not found in columns: {columns}"
        )
    if args.sort_column != "" and args.sort_column not in columns:
        raise ValueError(
            f"Sort column '{args.sort_column}' not found in columns: {columns}"
        )

    groups[first_row[args.query_column]].append(first_row)
    for row in rows_iter:
        groups[row[args.query_column]].append(row)

    # Sample
    talk_to_me(f"Sampling up to {args.n_alignments} alignments per query.")
    sampled_rows = []

    if args.sort_column == "":
        # Random mode
        random.seed(args.random_seed)
        for query_val, rows in groups.items():
            if len(rows) <= args.n_alignments:
                sampled_rows.extend(rows)
            else:
                sampled_rows.extend(random.sample(rows, args.n_alignments))
    else:
        # Top-N by sort_column (descending)
        for query_val, rows in groups.items():
            if len(rows) <= args.n_alignments:
                sampled_rows.extend(rows)
            else:
                sorted_rows = sorted(
                    rows,
                    key=lambda r: float(r[args.sort_column]),
                    reverse=True,
                )
                sampled_rows.extend(sorted_rows[: args.n_alignments])

    # Write output
    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)

    with open(args.output_file, "w") as outfile:
        outfile.write("\t".join(columns) + "\n")
        for row in sampled_rows:
            outfile.write("\t".join(row[c] for c in columns) + "\n")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
