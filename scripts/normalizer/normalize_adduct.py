import pandas as pd
import re
import os

global sub_adduct_pattern
sub_adduct_pattern = re.compile(r"\(|\)|(.*\[)|(\]([\d\+\-\*]*)?)")

global adduct_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_dict = dict(zip(adduct_dataframe['known_adduct'], adduct_dataframe['fraghub_default']))

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
    adduct = re.sub(sub_adduct_pattern, "", adduct)

    # Check if the normalized 'adduct' value exists in a previously defined dictionary 'adduct_dict'
    # If yes, replace the 'PRECURSORTYPE' value in the metadata dictionary with the corresponding value from 'adduct_dict'
    if adduct in adduct_dict:
        metadata_dict['PRECURSORTYPE'] = adduct_dict[adduct]

    # Return the metadata dictionary with normalized 'PRECURSORTYPE' value
    return metadata_dict
