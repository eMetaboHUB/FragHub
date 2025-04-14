import scripts.globals_vars
import re


def normalize_ms_level(metadata_dict):
    """
    Normalize the MS level value in the given metadata dictionary.
    :param metadata_dict: The dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized MS level value.
    :rtype: dict
    """
    try:
        ms_level = str(metadata_dict["MSLEVEL"])  # Retrieve MS level value from the given dictionary.
    except:
        ms_level = None

    if ms_level:  # Check if it exists or not None.
        ms_level = re.findall(globals_vars.ms_level_pattern, ms_level)  # Find all occurrences of the pattern in the MS level value.

        if ms_level:  # If pattern found,

            if len(ms_level) == 1:  # and if there's only one occurrence, set that occurrence as the MS level value.
                metadata_dict["MSLEVEL"] = ms_level[0]

            elif len(ms_level) >= 2:  # If there are two or more occurrences,
                # set the MS level value as a string that starts with the first occurrence and ends with the second one, separated by a dash.
                metadata_dict["MSLEVEL"] = f"{ms_level[0]}-{ms_level[1]}"

        else:  # If pattern not found, set MS level value as "2".
            metadata_dict["MSLEVEL"] = "2"  # Init MSLEVEL to 2 by default

    else:  # If MS level value does not exist or None, set its value as "2".
        metadata_dict["MSLEVEL"] = "2"  # Init MSLEVEL to 2 by default

    return metadata_dict  # Return the updated dictionary.
