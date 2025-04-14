import scripts.globals_vars
import numpy as np
import re


def normalize_empties(metadata_dict):
    """
    Normalizes empties in a metadata dictionary.
    :param metadata_dict: the dictionary containing metadata
    :return: the updated metadata dictionary
    """
    for k, v in metadata_dict.items():  # traversing all items (key-value pairs) in the dictionary

        if isinstance(v, str):  # if the value is a string
            if re.fullmatch(globals_vars.empty_pattern, v):  # if the value matches the 'empty_pattern' regex
                metadata_dict[k] = ''  # replace value in dictionary with empty string

        # if the value is a float or numpy float and is NaN,
        # replace the value in dictionary with empty string
        elif isinstance(v, (float, np.float64)) and np.isnan(v):
            metadata_dict[k] = ''

    return metadata_dict  # return updated dictionary
