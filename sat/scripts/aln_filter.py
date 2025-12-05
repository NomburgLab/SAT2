from .utils.misc import talk_to_me, make_output_dir


def validate_and_format_args(args):
    if args.alignment_fields != "":
        args.alignment_fields = args.alignment_fields.split(",")
    return args


def parse_filter_value(val_str):
    """Parse a string value to float, handling scientific notation."""
    try:
        return float(val_str.replace("E", "e"))
    except ValueError:
        raise ValueError(f"Cannot convert '{val_str}' to float for filtering.")


def aln_filter_main(args):
    args = validate_and_format_args(args)

    talk_to_me("Reading and filtering alignments.")

    alignment_fields = args.alignment_fields
    filter_field = args.filter_field
    min_val = args.min_val_filter_field
    max_val = args.max_val_filter_field

    if max_val < min_val:
        raise ValueError(
            f"max_val can't be less than min_val! max_val: {max_val}, min_val: {min_val}"
        )

    # List to store filtered alignments as (parts_list)
    filtered_alignments = []
    filter_field_idx = None

    with open(args.alignment_file) as infile:
        for line_num, line in enumerate(infile):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")

            # Handle header detection on first line
            if line_num == 0:
                if alignment_fields == "":
                    if parts[0] == "query":
                        alignment_fields = parts
                        continue
                    else:
                        raise ValueError(
                            "alignment_fields has not been provided, "
                            "which is only allowed when the first line has headers! "
                            "(e.g. first line should start with 'query')"
                        )
                else:
                    # If user provided fields, check if this is a header to skip
                    if parts[0] == "query":
                        continue

            # Skip header lines that might appear later in the file
            if parts[0] == "query":
                continue

            # Validate line length
            if len(parts) != len(alignment_fields):
                raise ValueError(
                    f"Line and alignment_fields don't have the same number of entries! "
                    f"Line has {len(parts)} entries, expected {len(alignment_fields)}. "
                    f"Line: {parts}"
                )

            # Find filter field index if not yet determined
            if filter_field_idx is None:
                try:
                    filter_field_idx = alignment_fields.index(filter_field)
                except ValueError:
                    raise ValueError(
                        f"Cannot find filter field '{filter_field}' in alignment fields: "
                        f"{alignment_fields}"
                    )

            # Get filter value and apply min/max filter
            filter_val = parse_filter_value(parts[filter_field_idx])

            if filter_val < min_val or filter_val > max_val:
                continue

            # Update the filter field with parsed float value
            parts[filter_field_idx] = str(filter_val)

            filtered_alignments.append(parts)

    talk_to_me("Writing output file.")
    make_output_dir(args.output_file)

    with open(args.output_file, "w") as outfile:
        # Write header
        outfile.write("\t".join(alignment_fields) + "\n")
        # Write filtered alignments
        for parts in filtered_alignments:
            outfile.write("\t".join(parts) + "\n")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
