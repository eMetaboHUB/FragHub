import re

global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)

def in_filename_or_name(filename, name):
    """
    Check whether a given filename or name matches specific criteria.

    This function checks if the string "MSMS_Public" is not in the filename
    and simultaneously if the filename combined with the name matches the
    In_Silico_pattern regular expression. If both conditions are met, the function
    returns True, otherwise, it returns False.

    Args:
        filename (str): The filename to be checked.
        name (str): The name to be checked alongside the filename.

    Returns:
        bool: True if conditions are met, False otherwise.
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
    Update the 'PREDICTED' field in the metadata dictionary based on specific conditions.

    The function checks if certain conditions are met within the metadata dictionary
    and updates the 'PREDICTED' field accordingly. The 'COMMENT' field is evaluated
    against a pattern, the 'PREDICTED' field is checked for a string "true", and the
    'FILENAME' field is checked for specific substrings.

    Parameters:
        metadata_dict (dict): A dictionary containing metadata with keys "COMMENT",
                              "PREDICTED", "FILENAME", and "NAME".

    Returns:
        dict: The updated metadata dictionary with the 'PREDICTED' field set to 'true'
              or 'false' based on the evaluated conditions.
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
