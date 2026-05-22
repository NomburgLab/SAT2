# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import os
from multiprocessing import Pool

from .utils.structure import pdb_to_structure_object, get_structure_paths, structure_to_pLDDT
from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def _compute_plddt(pdb_file):
    structure = pdb_to_structure_object(pdb_file)
    plddts = structure_to_pLDDT(structure)
    avg = round(sum(plddts.values()) / len(plddts), ndigits=2)
    return (os.path.basename(pdb_file), avg)




# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def struc_to_plddt_main(args):
    structure_paths = get_structure_paths(args.structure_file)
    if len(structure_paths) > 1:
        talk_to_me(f"Found {len(structure_paths)} PDB files in {args.structure_file}")

    threads = getattr(args, "threads", 1) or 1
    if threads > 1 and len(structure_paths) > 1:
        talk_to_me(f"Processing with {threads} workers")
        with Pool(processes=threads) as pool:
            results = pool.map(_compute_plddt, structure_paths)
    else:
        results = [_compute_plddt(f) for f in structure_paths]

    if args.out_file == "":
        for basename, plddt in results:
            if len(results) == 1:
                print(plddt)
            else:
                print(f"{basename}\t{plddt}")
    else:
        make_output_dir(args.out_file)
        with open(args.out_file, "a") as outfile:
            for basename, plddt in results:
                outfile.write(f"{basename}\t{plddt}\n")

    talk_to_me("Done!")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
