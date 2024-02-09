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

    :param metadata_dict: The dictionary containing metadata information.
    :return: The modified metadata dictionary with normalized adduct value.
    """
    adduct = metadata_dict['PRECURSORTYPE']
    adduct = re.sub(sub_adduct_pattern, "", adduct)

    if adduct in adduct_dict:
        metadata_dict['PRECURSORTYPE'] = adduct_dict[adduct]

    return metadata_dict