import pandas as pd
import scripts.globals_vars
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
    instrument_type = metadata_dict["INSTRUMENTTYPE"]
    if re.search("\b(GC)\b", instrument_type, flags=re.IGNORECASE):
        return metadata_dict

    # Get the 'PRECURSORTYPE' from the provided metadata dictionary
    adduct = metadata_dict['PRECURSORTYPE']

    # Normalize the obtained 'PRECURSORTYPE' value using regex substitution
    # Note: 'sub_adduct_pattern' is a previously defined regular expression pattern
    adduct = re.sub(globals_vars.sub_adduct_pattern, "", adduct)
    adduct = re.sub(globals_vars.sub_signe_end_adduct_pattern, "", adduct)

    # Check if the normalized 'adduct' value exists in a previously defined dictionary 'adduct_dict'
    # If yes, replace the 'PRECURSORTYPE' value in the metadata dictionary with the corresponding value from 'adduct_dict'
    if adduct in globals_vars.adduct_dict_POS:
        metadata_dict['PRECURSORTYPE'] = globals_vars.adduct_dict_POS[adduct]
    if adduct in globals_vars.adduct_dict_NEG:
        metadata_dict['PRECURSORTYPE'] = globals_vars.adduct_dict_NEG[adduct]

    # Return the metadata dictionary with normalized 'PRECURSORTYPE' value
    return metadata_dict
