import re

global ionization_mode_pattern
ionization_mode_pattern = re.compile(r"((?:^|\b)?APCI(?:\b|$)?)|((?:^|\b)?ACPI(?:\b|$)?)|((?:^|\b)?APPI(?:\b|$)?)|((?:^|\b)?EI(?:\b|$)?)|((?:^|\b)?ESI(?:\b|$)?)|((?:^|\b)?FAB(?:\b|$)?)|((?:^|\b)?MALDI(?:\b|$)?)",flags=re.IGNORECASE)

def normalize_ionization(metadata_dict):
    """
    Normalize the ionization mode in the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata.
    :type metadata_dict: dict
    :return: The modified metadata dictionary.
    :rtype: dict
    """
    ionization_mode = re.search(ionization_mode_pattern, metadata_dict["IONIZATION"])

    if ionization_mode:
        ionization_mode = ionization_mode.group(1)
        if ionization_mode == "ACPI":
            ionization_mode = "APCI"
        metadata_dict["IONIZATION"] = ionization_mode
    else:
        ionization_mode_in_INSTRUMENTTYPE = re.search(ionization_mode_pattern, metadata_dict["INSTRUMENTTYPE"])
        if ionization_mode_in_INSTRUMENTTYPE:
            ionization_mode_in_INSTRUMENTTYPE = ionization_mode_in_INSTRUMENTTYPE.group(1)
            if ionization_mode_in_INSTRUMENTTYPE == "ACPI":
                ionization_mode_in_INSTRUMENTTYPE = "APCI"
            metadata_dict["IONIZATION"] = ionization_mode_in_INSTRUMENTTYPE

    return metadata_dict