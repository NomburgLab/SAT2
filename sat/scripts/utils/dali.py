def parse_structure_key(structure_key_file, delim=",,", existing_dict=dict()):
    """
    Takes in a path to a file with the first column the structure, and second column
    the 4-digit identifier, and returns a dictionary of format
    identifier:structure.

    Add an existing_dict if you just want to update that dictionary - this lets you
    have multiple structure_key files.
    """

    key_to_structure = existing_dict.copy()
    with open(structure_key_file) as infile:
        for line in infile:
            line = line.rstrip("\n")
            structure, key = line.split(delim)

            if len(key) != 4:
                msg = "The key is expected to be a four-digit identifier! This key, "
                msg += f" {key}, is a length of {len(key)}!"
                raise ValueError(msg)

            if key in key_to_structure:
                msg = f"Have obseved a key, {key}, that is already present in"
                msg += " key_to_structure! This means it may be present multiple times!"
                raise ValueError(msg)

            key_to_structure[key] = structure

    return key_to_structure
