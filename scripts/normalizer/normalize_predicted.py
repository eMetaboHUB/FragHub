import re

global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)

def in_filename(filename):
    """
    :param filename: The name of the file to be checked for the presence of the string "MSMS_Public" and matching a specific pattern.
    :return: Returns True if the filename does not contain "MSMS_Public" and matches a specific pattern defined by the In_Silico_pattern. Returns False otherwise.

    """
    if "MSMS_Public" not in filename:
        if re.search(In_Silico_pattern, filename):
            return True

    return False

def normalize_predicted(metadata_dict):
    """
    Normalize the predicted field in the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata information.
    :return: The updated metadata dictionary with the normalized predicted field.
    """
    comment_field = metadata_dict["COMMENT"]
    predicted = metadata_dict["PREDICTED"]
    filename = metadata_dict["FILENAME"]

    if re.search(In_Silico_pattern, comment_field) or predicted == "true" or in_filename(filename):
        metadata_dict["PREDICTED"] = "true"
        return metadata_dict
    else:
        metadata_dict["PREDICTED"] = "false"
        return metadata_dict