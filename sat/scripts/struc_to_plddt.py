# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import os
from multiprocessing import Pool

from .utils.structure import pdb_to_structure_object, parse_structure_inputs, structure_to_pLDDT
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
    pdb_files = parse_structure_inputs(args.structure_file)
    if len(pdb_files) > 1:
        talk_to_me(f"Found {len(pdb_files)} PDB files in {args.structure_file}")

    threads = getattr(args, "threads", 1) or 1
    if threads > 1 and len(pdb_files) > 1:
        talk_to_me(f"Processing with {threads} workers")
        with Pool(processes=threads) as pool:
            results = pool.map(_compute_plddt, pdb_files)
    else:
        results = [_compute_plddt(f) for f in pdb_files]

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
