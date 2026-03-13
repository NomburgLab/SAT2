# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
from .utils.misc import talk_to_me, make_output_dir
import os

# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #

# Standard amino acid residue names found in PDB files
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    # Nonstandard but commonly accepted protein residues
    "MSE", "SEC", "PYL",
}

def ligands_cleaner(pdb_path):
    """
    Takes in a path to a PDB file and returns the lines cleaned of waters,
    ions, ligands, and any other non-protein residues. Only standard protein
    ATOM/TER/MODEL/ENDMDL/END records are kept; HETATM lines are retained
    only when the residue name is a recognised amino acid (e.g. MSE).
    """
    if not pdb_path.endswith(".pdb"):
        msg = "File must be a .pdb file and path must end with .pdb. pdb_path is "
        msg += f"currently set to {pdb_path}"
        raise ValueError(msg)
    if not os.path.exists(pdb_path):
        raise ValueError(f"Cannot detect PDB file: {pdb_path}")

    cleaned_lines = []
    with open(pdb_path) as infile:
        for line in infile:
            record = line[:6].strip()

            # Keep model / end bookkeeping lines as-is
            if record in ("MODEL", "ENDMDL", "END"):
                cleaned_lines.append(line)
                continue

            # For coordinate records, only keep standard amino acids
            if record in ("ATOM", "HETATM", "ANISOU"):
                resname = line[17:20].strip()
                if resname in STANDARD_AA:
                    cleaned_lines.append(line)
                continue

            # TER records: keep only if they follow a protein residue
            if record == "TER":
                # TER lines may or may not carry a residue name (cols 17-20)
                resname = line[17:20].strip() if len(line) > 20 else ""
                if resname == "" or resname in STANDARD_AA:
                    cleaned_lines.append(line)
                continue

            # Non-coordinate header/remark lines are dropped to give a
            # minimal, coordinates-only cleaned file.

    return cleaned_lines

# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
def struc_clean_main(args):
    talk_to_me("Reading PDB file")
    cleaned = ligands_cleaner(args.infile)

    if not cleaned:
        raise RuntimeError(
            f"No protein residues found in {args.infile}. The cleaned output "
            "would be empty — check that the input is a valid PDB with protein chains."
        )

    talk_to_me("Writing output.")
    make_output_dir(args.outfile)

    if os.path.exists(args.outfile):
        msg = "Output file already exists and will be overwritten: "
        msg += f"{args.outfile}"
        talk_to_me(msg)

    with open(args.outfile, "w") as outfile:
        outfile.writelines(cleaned)

    talk_to_me(f"Done. Cleaned PDB written to {args.outfile}")


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)