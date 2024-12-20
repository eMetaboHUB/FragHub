import globals_vars
import re


def normalize_ionization(metadata_dict):
    """
    Normalize the ionization mode in the given metadata dictionary.
    This function conducts two checks in order to normalize the 'IONIZATION' in the metadata.
    The first check is performed on 'IONIZATION' directly.
    If 'IONIZATION' is found, it's processed and the processed value is placed back into 'IONIZATION'.
    If 'IONIZATION' is not fpound, then 'INSTRUMENTTYPE' key is checked for ionization mode.
    If ionization mode is found in 'INSTRUMENTTYPE', it's then re-assigned to 'IONIZATION' key.
    If ionization mode is not found in any of the keys, 'IONIZATION' is set to empty string.
    Finally, the modified dictionary is returned.

    :param metadata_dict: A dictionary containing metadata.
    :type metadata_dict: dict
    :return: The modified metadata dictionary with normalized ionization mode.
    :rtype: dict
    """
    ionization_mode = re.search(globals_vars.ionization_mode_pattern, metadata_dict["IONIZATION"])
    if ionization_mode:
        ionization_mode = ionization_mode.group(1)
        if ionization_mode == "ACPI":  # correct known typo in ionization mode spelling
            ionization_mode = "APCI"
        metadata_dict["IONIZATION"] = ionization_mode
    else:
        ionization_mode_in_INSTRUMENTTYPE = re.search(globals_vars.ionization_mode_pattern, metadata_dict["INSTRUMENTTYPE"])
        if ionization_mode_in_INSTRUMENTTYPE:
            ionization_mode_in_INSTRUMENTTYPE = ionization_mode_in_INSTRUMENTTYPE.group(1)
            if ionization_mode_in_INSTRUMENTTYPE == "ACPI":  # correct known typo in ionization mode spelling
                ionization_mode_in_INSTRUMENTTYPE = "APCI"
            metadata_dict["IONIZATION"] = ionization_mode_in_INSTRUMENTTYPE
    if metadata_dict["IONIZATION"] == None:  # if no ionization mode is found, set the ionization field to empty string
        metadata_dict["IONIZATION"] = ''
    return metadata_dict
    # The modified dictionary is returned with normalized ionization mode
