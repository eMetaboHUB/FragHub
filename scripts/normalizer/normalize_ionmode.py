import re

global ionmode_pos_pattern
ionmode_pos_pattern = re.compile(r"^p|^\+|^pos", flags=re.IGNORECASE)

global ionmode_neg_pattern
ionmode_neg_pattern = re.compile(r"^n|^\-|^neg", flags=re.IGNORECASE)

def normalize_ionmode(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
        It should have the key "IONMODE" to represent the ionization mode.
    :return: The updated metadata dictionary with the "IONMODE" value normalized to either "positive" or "negative".
    """
    ionmode = metadata_dict["IONMODE"]

    if re.search(ionmode_pos_pattern, ionmode):
        ionmode = "positive"
    elif re.search(ionmode_neg_pattern, ionmode):
        ionmode = "negative"

    metadata_dict["IONMODE"] = ionmode

    return metadata_dict