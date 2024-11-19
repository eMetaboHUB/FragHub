import re

global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)

def in_filename_or_name(filename, name):
    """
    This function checks if a provided filename is valid based on two criteria:
    1. The filename does not contain the string "MSMS_Public"
    2. The filename matches a specific pattern (In_Silico_pattern) defined elsewhere in the code.

    :param filename: The name of the file to be checked for the presence of the string "MSMS_Public" and matching a specific pattern.
    :return: The function returns True if the filename meets both criteria. If the filename contains "MSMS_Public" or does not match the In_Silico_pattern, the function returns False.
    """
    # Check if the string "MSMS_Public" is not in the filename
    if "MSMS_Public" not in filename:
        # If "MSMS_Public" is not in the filename, check if it matches the In_Silico_pattern
        if re.search(In_Silico_pattern, filename+" "+name):
            # If both conditions are met, return True
            return True
    # If either of the conditions is not met, return False
    return False

def normalize_predicted(metadata_dict):
    """
    Normalize the predicted field in the given metadata dictionary.

    This function checks the 'COMMENT' and 'PREDICTED' values in the metadata dictionary and the filename
    of the data. If the 'COMMENT' field matches a certain pattern (represented by a regular expression
    'In_Silico_pattern'), or the 'PREDICTED' value is 'true', or the string "MSMS_Public" is in the filename,
    it sets the 'PREDICTED' value to 'true', else it sets the 'PREDICTED' value to 'false'.

    :param metadata_dict: A dictionary containing metadata information.
    :return: The updated metadata dictionary with the normalized predicted field.
    """

    comment_field = metadata_dict["COMMENT"]  # Extract the 'COMMENT' field from the metadata dictionary
    predicted = metadata_dict["PREDICTED"]  # Extract the 'PREDICTED' field from the metadata dictionary
    filename = metadata_dict["FILENAME"]  # Extract the 'FILENAME' from the metadata dictionary
    name = metadata_dict["NAME"]  # Extract the 'NAME' from the metadata dictionary

    # If 'COMMENT' field matches the pattern, or 'PREDICTED' is 'true', or 'MSMS_Public' in the filename:
    #    set 'PREDICTED' field in the metadata dictionary to 'true'
    if re.search(In_Silico_pattern, comment_field) or predicted == "true" or in_filename_or_name(filename, name):
        metadata_dict["PREDICTED"] = "true"
    else:  # Otherwise, set the 'PREDICTED' field in the metadata dictionary to 'false'
        metadata_dict["PREDICTED"] = "false"

    return metadata_dict  # Return the updated metadata dictionary
