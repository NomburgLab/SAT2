import os
import json

from .utils.misc import talk_to_me, make_output_dir


def get_iptm_from_json(json_path):
    """
    Takes in path to a colabfold json file, and takes the iptm out of there.

    Note that even if AF2 is run with the monomer model, as long as the input was
    a multimer there will still be an iptm in the json.
    """
    if not json_path.endswith(".json"):
        msg = "File must be a .json file and path must end with .json. json_path is "
        msg += f"currently set to {json_path}"
        raise ValueError(msg)

    if not os.path.exists(json_path):
        raise ValueError(f"Cannot detect json file: {json_path}")

    with open(json_path) as infile:
        data = json.load(infile)
        if "iptm" not in data:
            msg = f"Cannot find iptm in the input json, {json_path}. It is expected "
            msg += "from both AF2 monomer and multimer."
            raise ValueError(msg)
        return str(data["iptm"])


def struc_get_iptm_main(args):
    talk_to_me("Reading json file")
    iptm = get_iptm_from_json(args.infile)
    file_basename = os.path.basename(args.infile.rstrip(".json"))
    out_line = f"{file_basename}\t{iptm}\n"

    talk_to_me("Writing output.")
    make_output_dir(args.outfile)
    if os.path.exists(args.outfile):
        msg = "This script appends to the output file. An output file already exists, so it"
        msg += " will be appended to!"
        talk_to_me(msg)
    with open(args.outfile, "a") as outfile:
        outfile.write(out_line)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
