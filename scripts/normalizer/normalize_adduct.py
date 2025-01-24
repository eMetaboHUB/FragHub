import deletion_report
import pandas as pd
import globals_vars
import re
import os


def normalize_adduct(metadata_dict):
    """
    Normalize adduct value in the given metadata dictionary.
    The function is getting the 'PRECURSORTYPE' from metadata dictionary,
    normalize it and replacing it back in the metadata dictionary.
    :param metadata_dict: The dictionary containing metadata information.
    :return: The modified metadata dictionary with normalized adduct value.
    """
    # Get the 'PRECURSORTYPE' from the provided metadata dictionary
    adduct = metadata_dict['PRECURSORTYPE']

    # Normalize the obtained 'PRECURSORTYPE' value using regex substitution
    # Note: 'sub_adduct_pattern' is a previously defined regular expression pattern
    adduct = re.sub(globals_vars.sub_adduct_pattern, "", adduct)

    if not re.search(globals_vars.is_adduct_pattern, adduct):
        deletion_report.no_or_bad_adduct += 1
        return None

    # Check if the normalized 'adduct' value exists in a previously defined dictionary 'adduct_dict'
    # If yes, replace the 'PRECURSORTYPE' value in the metadata dictionary with the corresponding value from 'adduct_dict'
    if adduct in globals_vars.adduct_dict:
        metadata_dict['PRECURSORTYPE'] = globals_vars.adduct_dict[adduct]

    # Return the metadata dictionary with normalized 'PRECURSORTYPE' value
    return metadata_dict
