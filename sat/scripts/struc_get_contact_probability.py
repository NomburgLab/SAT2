# ------------------------------------------------------------------------------------ #
# Import dependencies
# ------------------------------------------------------------------------------------ #
import pickle
import numpy as np
import scipy.special
from glob import glob
import os
import re
import json

from .utils.misc import talk_to_me, make_output_dir


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def read_distogram(pkl_file):
    """
    Input: path to pickle file
    Output: Distogram
    """
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    distogram = data.get("distogram", None)
    if distogram:
        return distogram
    else:
        msg = "Cannot detect a distogram within the pickle object. Something is wrong."
        raise ValueError(msg)


def sum_contact_probabilities(distogram, length_cutoff):
    """
    Summarize contact probabilities for residue pairs within a given distance cutoff.

    Parameters:
    - distogram (dict): AlphaFold distogram object containing "logits" and "bin_edges".
    - length_cutoff (float): Distance cutoff (e.g., 8.0 or 12.0 Å).

    Returns:
    - np.ndarray: A 2D array (N x N) of summed contact probabilities for residue pairs.

    More lay-man description. This function does the following:
    - The input is the distogram, which has an array of bin edges and a 3D matrix of
      logits.

      The bin edges correspond to different distance cutoffs bins. There are D bins
      with a cutoff (e.g. 2A, 6A, 9A,...)

      The logits matrix has shape N x N x D + 1, where N is equal to the total
      number of residues in the input structure, while D is equal to the total number of
      distance bins (excluding the outlier distance bin). Basically, NxN array 1
      corresponds to logits from the smallest distance bin, NxN array 2 is the next one,
      etc etc. The extra dimension of D is the logits (that could be converted to
      probability) of the residue-residue contact being above the final distance bin.

    - This function first uses the softmax function to convert logits to probability.
      Basically, each cell is exponentialized with exp() (solve for e^value). Then,
      each of those values is normalized across all of the D+1 bins - e.g. the first
      row/first col cell of the first bin is normalized to the first row/first col
      cell of all other bins. Now, you still have a NxNxD matrix, but each value is the
      probability that the residue-residue contact is within that distance bin.

    - This function goes a final step - it sums all of the probabilities for each
      residue-residue cell for all bins within the distance cutoff. So, you first find
      using the bin edges which bins are below distance cutoff, then sum those. This
      in essence is finding the probability for each residue-residue cell that it
      falls below that distance cutoff. So, the final output is an NxN matrix.
    """
    # Extract logits and bin_edges from distogram
    logits = distogram.get("logits")
    bin_edges = distogram.get("bin_edges")

    if logits is None or bin_edges is None:
        raise ValueError("Distogram must contain 'logits' and 'bin_edges'.")

    # Convert logits to probabilities using softmax
    probabilities = scipy.special.softmax(logits, axis=-1)  # Shape: (N, N, D+1)

    # Identify bins that correspond to distances ≤ length_cutoff
    contact_bins = np.where(bin_edges <= length_cutoff)[0]

    if len(contact_bins) == 0:
        raise ValueError(f"No bins found for the given length cutoff: {length_cutoff}")

    # Sum probabilities across the selected bins
    summed_probabilities = np.sum(
        probabilities[:, :, contact_bins], axis=-1
    )  # Shape: (N, N)

    return summed_probabilities


def highest_contact_probability(summed_probabilities, len1):
    """
    Extract the submatrix representing contacts between Protein 1 and Protein 2 and return the highest contact probability.

    Parameters:
    - summed_probabilities (np.ndarray): A 2D array (N x N) of contact probabilities.
    - len1 (int): Length of Protein 1. Protein 2 is assumed to have length (N - len1).

    Returns:
    - float: The highest contact probability in the submatrix.

    More lay-man description. This function does the following:
    - previously, we converted logits to probabilities, then summed the probabilities
      for bins within the indicated distance cutoff. So, the input summed_probabilities
      is an NxN matrix. Here, we are finding the submatrix m' that corresponds to
      the residue-residue contacts between protein chains. This evaluates to
      m' = m[: len1][len1:len1 + len2]. This is bascially all of the residue-residue
      contacts between chains.
    - this function extracts that submatrix, and then finds the highest
      value in that submatrix. This is the highest residue-residue contact probability
      between the two chains.
    - See Figure S6 of 10.1126/science.abm4805 for a visual representation/source
      methods.
    """
    # Total length of the summed_probabilities matrix
    total_len = summed_probabilities.shape[0]

    # Ensure len1 is valid
    if len1 <= 0 or len1 >= total_len:
        raise ValueError(
            "Invalid Protein 1 length (len1). Ensure 0 < len1 < total_len."
        )

    # Define the length of Protein 2
    len2 = total_len - len1

    # Extract the submatrix m'
    submatrix = summed_probabilities[:len1, len1 : len1 + len2]

    # Find the highest contact probability in the submatrix
    max_contact_probability = np.max(submatrix)

    return max_contact_probability


def parse_pr1_len(a3m_file_path):
    """
    This function parses the colabfold a3m file to retreive protein 1 length. This
    is present as the first line of the colabfold a3m file, which looks like:
    #98,76  1,1
    Where this input file contained two proteins, one 98 residues and one 76 residues.

    This function would return 98 in this example.
    """

    with open(a3m_file_path) as infile:
        first_line = infile.readline().strip()
        pr1_len = int(first_line.split("\t")[0].lstrip("#").split(",")[0])
        return int(pr1_len)


def pickle_name_to_info(input):
    """
    Takes in the name of a pickle file, and parses the prefix, model type, rank, and model.
    For example
    input: 6A6I_1__6H4B_1_all_rank_001_alphafold2_ptm_model_4_seed_000
    output: A tuple of 6A6I_1__6H4B_1, alphafold2_ptm, 1, 4
    """
    if "_all_rank" not in input:
        msg = (
            "Cannot find the substring '_all_rank' in the name of the pickle, "
            f"{input}. Something is wrong."
        )
        raise ValueError(msg)
    sample_name = input.split("_all_rank")[0]
    rank_match = re.search(r"all_rank_(\d+)", input)
    model_match = re.search(r"model_(\d+)", input)
    model_type_match = re.search(r"all_rank_\d+_(.*?)_model_\d+", input)

    if not rank_match or not model_match or not model_type_match:
        msg = "Cannot parse the strings all_rank_ and/or model_ from the pickle file name, "
        msg += f"{input}. Something is wrong."
        raise ValueError(msg)

    rank = rank_match.group(1).lstrip("0")
    model = model_match.group(1)
    model_type = model_type_match.group(1)

    return sample_name, model_type, rank, model


def get_iptm_from_json(json_path):
    """
    Takes in path to a colabfold json file, and takes the iptm out of there.

    iPTM should be present for both AF monomer and mulitmer.
    """
    if not os.path.exists(json_path):
        raise ValueError(f"Cannot detect json file: {json_path}")

    with open(json_path) as infile:
        data = json.load(infile)
        if "iptm" not in data:
            msg = f"Cannot find iptm in the input json, {json_path}. It is expected "
            msg += "from both AF2 monomer and multimer."
            raise ValueError(msg)
        return str(data["iptm"])


def struc_get_contact_probability_main(args):
    talk_to_me("Getting protein length from a3m")
    a3m = glob(f"{args.in_dir}/*a3m")
    if len(a3m) > 1:
        raise ValueError("There are more than one a3m files in the directory, weird.")
    pr1_len = parse_pr1_len(a3m[0])

    talk_to_me("Iterating over pickle files")
    output = ""
    for FILE in glob(f"{args.in_dir}/*pickle"):
        basename = os.path.basename(FILE).replace(".pickle", "")
        sample_name, model_type, rank, model = pickle_name_to_info(basename)

        # Parsring the distogram from the pickle
        distogram = read_distogram(FILE)
        contact_prbabilities = sum_contact_probabilities(
            distogram, args.distance_cutoff
        )
        prob = str(highest_contact_probability(contact_prbabilities, pr1_len))

        # Also parse the iPTM from the json file corresponding to each pickle
        json_path = f"{FILE.rstrip('.pickle')}.json".replace("all", "scores")
        iptm = get_iptm_from_json(json_path)

        out_line = "\t".join([sample_name, model_type, rank, model, iptm, prob]) + "\n"
        output += out_line

    talk_to_me("Writing output")
    make_output_dir(args.outfile)
    with open(args.outfile, "w") as outfile:
        outfile.write(output)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
