import re

def normalize_retentiontime(metadata_dict):
    """
    Normalize the retention time value in the given metadata dictionary.
    :param metadata_dict: The dictionary containing the metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized retention time.
    :rtype: dict
    """

    # Retrieving the retention time from metadata dictionary.
    try:
        retientiontime = str(metadata_dict["RETENTIONTIME"])
    except:
        retientiontime = ""

    # Regular expression match for different patterns of time representation.
    match = re.search(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(m|min|minute|minutes|s|sec|second|seconds|ms|millisecond|milliseconds)(?:\W)?", retientiontime, flags=re.IGNORECASE)

    # If match found, normalization process begins.
    if match:
        # Group 1 in regex match refers to actual time representation.
        time = match.group(1)

        # Group 2 in regex match refers to the units (minutes, seconds, milliseconds).
        unit = match.group(2).lower()

        # If no unit specification found, default unit is considered as minutes (m).
        # The time value is converted to float and stored as string in 'RETENTIONTIME'.
        if not unit:
            retientiontime = str(float(time))
            metadata_dict["RETENTIONTIME"] = retientiontime
            return metadata_dict
        else:
            # If the unit is minutes.
            # Example: 'm', 'min', 'minute', 'minutes', then directly store as string, considering no normalization is required.
            if unit in ["m", "min", "minute", "minutes"]:
                retientiontime = str(float(time))
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict

            # If the unit is seconds. Example: 's', 'sec', 'second', or 'seconds', then convert it to minutes by dividing by 60.
            elif unit in ["s", "sec", "second", "seconds"]:
                retientiontime = str(float(time) / 60)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict

            # If the unit is milliseconds. Example: 'ms', 'millisecond', or 'milliseconds', then convert it to minutes by dividing 60000.
            elif unit in ["ms", "millisecond", "milliseconds"]:
                retientiontime = str(float(time) / 60000)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict

    # If no match found in the above cases, return the metadata dictionary as it is.
    return metadata_dict
