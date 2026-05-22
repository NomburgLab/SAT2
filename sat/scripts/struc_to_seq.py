# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import os
from pathlib import Path
from multiprocessing import Pool

from .utils.structure import pdb_to_structure_object, get_structure_paths, struc_to_seq
from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def _extract_seq(pdb_file):
    basename = os.path.basename(pdb_file).removesuffix(".pdb")
    try:
        structure = pdb_to_structure_object(pdb_file, basename)
        seq = struc_to_seq(structure)
    except Exception as e:
        raise ValueError(f"Failed to parse PDB file: {pdb_file}") from e
    if not seq:
        raise ValueError(f"Failed to parse PDB file: {pdb_file}")
    return f">{basename}\n{seq}\n"



# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def struc_to_seq_main(args):
    structure_paths = get_structure_paths(args.structure_file)
    if len(structure_paths) > 1:
        talk_to_me(f"Found {len(structure_paths)} PDB files in {args.structure_file}")
    directory_mode = len(structure_paths) > 1 or Path(args.structure_file).is_dir()

    if directory_mode:
        threads = getattr(args, "threads", 1) or 1
        if threads > 1:
            talk_to_me(f"Processing with {threads} workers")
            with Pool(processes=threads) as pool:
                results = pool.imap(_extract_seq, structure_paths, chunksize=50)
                _write_directory_results(results, len(structure_paths), args)
        else:
            results = (_extract_seq(f) for f in structure_paths)
            _write_directory_results(results, len(structure_paths), args)
    else:
        _write_single_file(structure_paths[0], args)

    talk_to_me("Done!")


def _write_directory_results(results, total, args):
    n_ok = 0

    if args.out_file == "":
        for i, seq_record in enumerate(results):
            print(seq_record, end="")
            n_ok += 1
            if (i + 1) % 10000 == 0:
                talk_to_me(f"  {i + 1}/{total} done")
    else:
        make_output_dir(args.out_file)
        with open(args.out_file, "w") as out_f:
            for i, seq_record in enumerate(results):
                out_f.write(seq_record)
                n_ok += 1
                if (i + 1) % 10000 == 0:
                    talk_to_me(f"  {i + 1}/{total} done")

    talk_to_me(f"Wrote {n_ok} sequences")


def _write_single_file(pdb_file, args):
    header = getattr(args, "header", "") or ""
    structure = pdb_to_structure_object(pdb_file, header or "structure")
    seq = struc_to_seq(structure)

    if args.out_file == "":
        print(seq)
    else:
        if header == "":
            msg = "If --out_file is specified, you must fill in --header."
            raise ValueError(msg)
        make_output_dir(args.out_file)
        with open(args.out_file, "a") as outfile:
            out = f">{header}\n{seq}\n"
            outfile.write(out)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
