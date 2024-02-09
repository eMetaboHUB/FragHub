import re

global ms_level_pattern
ms_level_pattern = re.compile(r"(?:ms)?(\d)", flags=re.IGNORECASE)

def normalize_ms_level(metadata_dict):
    """
    Normalize the MS level value in the given metadata dictionary.

    :param metadata_dict: The dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized MS level value.
    :rtype: dict
    """
    ms_level = metadata_dict["MSLEVEL"]
    if ms_level:
        ms_level = re.findall(ms_level_pattern, ms_level)
        if ms_level:
            if len(ms_level) == 1:
                metadata_dict["MSLEVEL"] = ms_level[0]
            elif len(ms_level) >= 2:
                metadata_dict["MSLEVEL"] = f"{ms_level[0]}-{ms_level[1]}"
        else:
            metadata_dict["MSLEVEL"] = "2" # Init MSLEVEL to 2 by default
    else:
        metadata_dict["MSLEVEL"] = "2"  # Init MSLEVEL to 2 by default

    return metadata_dict