import numpy as np
import re

def normalize_empties(metadata_dict):
    regex_to_replace = [re.compile(re.escape(str(item)), re.I) for item in [0, 0.0, "", "na", "n/a", "nan", "unknown","unknown","none", np.nan]]

    for k, v in metadata_dict.items():
        # For each item to replace
        for regex in regex_to_replace:
            # If the current value matches, replace it
            if regex.fullmatch(str(v)):
                metadata_dict[k] = np.nan

    return metadata_dict

def normalize_values(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    metadata_dict = normalize_empties(metadata_dict)
    
    # metadata_dict = normalize_adduct(metadata_dict)
    # metadata_dict = normalize_ionmode(metadata_dict)
    # metadata_dict = normalize_retention_time(metadata_dict)
    # metadata_dict = normalize_ms_level(metadata_dict)
    # metadata_dict = normalize_synonymes(metadata_dict)
    # metadata_dict = normalize_formula(metadata_dict)
    # metadata_dict = normalize_predicted(metadata_dict)
    # metadata_dict = normalize_db_informations(metadata_dict)

    return metadata_dict