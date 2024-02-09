import re

def normalize_retentiontime(metadata_dict):
    """
    Normalize the retention time value in the given metadata dictionary.

    :param metadata_dict: The dictionary containing the metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized retention time.
    :rtype: dict
    """
    retientiontime = metadata_dict["RETENTIONTIME"]

    match = re.search(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(m|min|minute|minutes|s|sec|second|seconds|ms|millisecond|milliseconds)(?:\W)?", retientiontime, flags=re.IGNORECASE)

    if match:
        time = match.group(1)
        unit = match.group(2).lower()

        if not unit:
            retientiontime = str(float(time))
            metadata_dict["RETENTIONTIME"] = retientiontime
            return metadata_dict
        else:
            if unit in ["m", "min", "minute", "minutes"]:
                retientiontime = str(float(time))
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict
            elif unit in ["s", "sec", "second", "seconds"]:
                retientiontime = str(float(time) / 60)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict
            elif unit in ["ms", "millisecond", "milliseconds"]:
                retientiontime = str(float(time) / 60000)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict

    return metadata_dict