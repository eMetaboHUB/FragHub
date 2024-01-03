import re

metadata_dict = {}
metadata_dict["RETENTIONTIME"] = "10.517 min"

def normalize_retentiontime(metadata_dict):
    retientiontime = metadata_dict["RETENTIONTIME"]

    match = re.search("(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(.*)", retientiontime, flags=re.IGNORECASE)

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

print(normalize_retentiontime(metadata_dict))
