import globals_vars
import json
import re


def clean_instrument(instrument):
    """
    Clean the instrument name string by removing specific prefixes and modifying the format.
    :param instrument: The instrument name string to be cleaned.
    :return: The cleaned instrument name string.
    """
    # Replace "-tof" with "tof" in instrument name
    instrument = re.sub("-tof", "tof", instrument)

    # Replace "q-" with "q" in instrument name
    instrument = re.sub("q-", "q", instrument)

    # Replace "q exactive" with " qexactive " in instrument name
    instrument = re.sub("q exactive", " qexactive ", instrument)

    # Replace "applied biosystems" with " sciex " in instrument name
    instrument = re.sub("applied biosystems", " sciex ", instrument)

    # Replace " ab " with " sciex " in instrument name
    instrument = re.sub(" ab ", " sciex ", instrument)

    # Add spaces around "sciex" in instrument name
    instrument = re.sub("sciex", " sciex ", instrument)

    # Replace "triple(-| )?tof" or "triple(-| )?tof" with " qqq " in instrument name
    instrument = re.sub("triple(-| )?tof", " qqq ", instrument)
    instrument = re.sub("triple(-| )?quad", " qqq ", instrument)

    # Remove "... uplc ..." from instrument name
    instrument = re.sub("... uplc ...", " ", instrument)

    # Return the cleaned instrument name
    return instrument

def clean_instrument_type(instrument_type):
    """
    Cleans the given instrument string by removing specific substrings and replacing hyphens with spaces.
    :param instrument_str: The instrument string to be cleaned.
    :type instrument_str: str
    :return: The cleaned instrument string.
    :rtype: str
    """
    print("before clean_instrument_type: ", instrument_type)
    # If the given instrument type string contains "-tof", this line replaces
    # it with "tof"
    instrument_type = re.sub("-tof", "tof", instrument_type)

    # If the given instrument type string contains "q-", this line replaces
    # it with "q"
    instrument_type = re.sub("q-", "q", instrument_type)

    # Replaces all hyphen characters in the instrument type string with spaces
    instrument_type = re.sub("-", " ", instrument_type)

    # If the given instrument type string contains "q exactive", this line
    # replaces it with "qexactive"
    instrument_type = re.sub("q exactive", " qexactive ", instrument_type)

    # If the given instrument type string contains "applied biosystems",
    # this line replaces it with "sciex"
    instrument_type = re.sub("applied biosystems", " sciex ", instrument_type)

    # If the given instrument type string contains " ab ", this line replaces
    # it with "sciex"
    instrument_type = re.sub(" ab ", " sciex ", instrument_type)

    # Surrounds "sciex" in the instrument type string with spaces, if not surrounded yet.
    instrument_type = re.sub("sciex", " sciex ", instrument_type)

    # If the given instrument type string contains "triple(-| )?tof" or
    # "triple(-| )?quad", this line replaces it with "qqq"
    instrument_type = re.sub("triple(-| )?tof", " qqq ", instrument_type)
    instrument_type = re.sub("triple(-| )?quad", " qqq ", instrument_type)

    # If the given instrument type string contains "... uplc ...", this line removes it
    instrument_type = re.sub("... uplc ...", " ", instrument_type)
    print("after clean_instrument: ", instrument_type)

    return instrument_type

def clean_comment(comment):
    """
    Cleans the given comment by removing specific strings and replacing certain characters.
    :param comment: The comment to be cleaned.
    :type comment: str
    :return: The cleaned comment.
    :rtype: str
    """

    # Replaces "-tof" with "tof" in the comment string
    comment = re.sub("-tof", "tof", comment)

    # Replaces "q-" with "q" in the comment string
    comment = re.sub("q-", "q", comment)

    # Replaces "-" (hyphen) with " " (space) in the comment string
    comment = re.sub("-", " ", comment)

    # Replaces "q exactive" with " qexactive " in the comment string
    comment = re.sub("q exactive", " qexactive ", comment)

    # Replaces "applied biosystems" with " sciex " in the comment string
    comment = re.sub("applied biosystems", " sciex ", comment)

    # Replaces " ab " with " sciex " in the comment string
    comment = re.sub(" ab ", " sciex ", comment)

    # Surrounds "sciex" in the comment string with spaces, if not surrounded yet
    comment = re.sub("sciex", " sciex ", comment)

    # Replaces "triple(-| )?tof" or "triple(-| )?quad" with " qqq " in the comment string
    comment = re.sub("triple(-| )?tof", " qqq ", comment)
    comment = re.sub("triple(-| )?quad", " qqq ", comment)

    # Removes "... uplc ..." from the comment string
    comment = re.sub("... uplc ...", " ", comment)

    # Returns the cleaned comment string
    return comment

def clean_spectrum_instrument_info(metadata_dict):
    """
    Cleans the spectrum instrument information from the given metadata dictionary.
    :param metadata_dict: The dictionary containing the metadata information.
    :return: The cleaned instrument information.
    """
    # extract instrument information from the metadata and convert to lowercase
    instrument = metadata_dict['INSTRUMENT'].lower()
    instrument_type = metadata_dict["INSTRUMENTTYPE"].lower()
    comment = metadata_dict["COMMENT"].lower()

    # clean individual instrument information
    instrument = clean_instrument(instrument)
    instrument_type = clean_instrument_type(instrument_type)
    comment = clean_comment(comment)

    # combine the cleaned instrument info
    instrument_infos = instrument + " " + instrument_type + " " + comment

    # remove any non-word, non-whitespace, and non-hyphen characters, and strip leading/trailing whitespaces
    instrument_infos = re.sub(r'[^-\w\s]', ' ', instrument_infos)
    instrument_infos = ' '.join(instrument_infos.split()).strip()

    return instrument_infos

def search_for_brand(tree_path, instrument_infos):
    """
    Searches for a brand in instrument_infos and appends the found brand or 'unknown' to the tree_path.
    :param tree_path: A list representing the path to the brand in the instrument_tree.
    :param instrument_infos: A string containing information about the instrument.
    :return: The updated tree_path list.
    """

    # try to execute the following instruction
    try:
        # For each key in the instrument tree dictionary
        for key in globals_vars.instrument_tree.keys():
            # If the key (brand) is found in the instrument_infos string
            if re.search(rf"(\b|^|$){key}(\b|^|$)", instrument_infos):
                # Append the found key (brand) to the tree_path list
                tree_path.append(key)
                # Return the updated tree_path list
                return tree_path
        # If no brand key is found, append 'unknown' to the tree_path list
        tree_path.append('unknown')
        # Return the updated tree_path list with 'unknown' appended
        return tree_path
    # If the try block throws an exception, return None
    except:
        return None

def search_for_model(tree_path, instrument_infos):
    """
    Search for a model in the instrument tree based on the given tree path and instrument infos.
    :param tree_path: The path to the model in the instrument tree.
    :param instrument_infos: The instrument infos used to search for the model.
    :return: The updated tree path with the found model or 'unknown' if not found.
    """

    # Start by assuming no exception will occur
    try:
        # Loop through the keys present in the first level of instrument_tree dictionary
        for key in globals_vars.instrument_tree[tree_path[0]].keys():
            # Using a regular expression, search for the key in instrument_infos.
            # rf"(\b){key}(\b)" is a pattern that matches the key surrounded by word boundaries to avoid partial matches.
            if re.search(rf"(\b|^|$){key}(\b|^|$)", instrument_infos):
                # If a match is found, append the key to the tree_path list
                tree_path.append(key)
                # Return the updated tree_path
                return tree_path

        # If no match is found after the loop, append 'unknown' to the tree_path list
        tree_path.append('unknown')

        # Return the updated tree_path
        return tree_path
    # If an exception occurs (like a KeyError), handle it by returning None
    except:
        return None

def search_for_spectrum_type(tree_path, instrument_infos):
    """
    Method to search for spectrum type in a given tree path and instrument infos.
    :param tree_path: A List of keys that represent a path in a nested dictionary structure.
    :param instrument_infos: A String which contains relevant instrument information.

    This function tries to find a matching key in the dictionary structure at the first two levels of the tree path
    and if found that key is append to the tree_path and returned.
    It does this by iterating over all key in that level of the dictionary structure and uses
    a regular expression to check if the key exists in the instrument information string.
    If no key is found the function appends 'unknown' to the tree_path and returns it.
    Should anything unexpected occur the function will return None.

    :return: The tree_path List with found spectrum type appended, or 'unknown' if not found. None if an exception occurred.
    """
    try:  # Try to execute the following block of code
        for key in globals_vars.instrument_tree[tree_path[0]][tree_path[1]].keys():  # This block of code iterates over each key found in the 2nd depth of the dictionary
            if re.search(rf"(\b|^|$){key}(\b|^|$)", instrument_infos):  # If the key is found in the instrument data string
                tree_path.append(key)  # Then append the found key value to the tree path
                return tree_path  # Return the updated tree path
        tree_path.append('unknown')  # If the key is not found in the instrument data string, append 'unknown' to the tree path
        return tree_path  # Return the updated tree path
    except:  # If any error occur while executing the above block of code
        return None  # Return None

def search_for_instrument_type(tree_path, instrument_infos):
    """
    Search for instrument type in instrument_infos based on tree_path.
    :param tree_path: list of indices representing the path in the instrument tree
    :type tree_path: list[int]
    :param instrument_infos: information about instruments
    :type instrument_infos: str
    :return: updated tree_path with the identified instrument type appended, otherwise 'unknown' if not found
    :rtype: list[int]
    """
    # Trying to match and append an element of the path to the tree_path,
    # If there is any error occurs, it will return None
    try:

        # Loop through the keys in the instrument tree at the provided tree path
        for key in globals_vars.instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]].keys():

            # Use regex to find if the current key exists in instrument_infos
            # If it does, append it to the tree_path and return the updated tree_path
            if re.search(rf"(\b|^|$){key}(\b|^|$)", instrument_infos):
                tree_path.append(key)
                return tree_path

        # If we've looped through every key and haven't returned,
        # append 'unknown' to the tree_path and return the updated tree_path
        tree_path.append('unknown')
        return tree_path

    except:
        # If any exception occurs, return None
        return None

def search_for_ionisation(tree_path, instrument_infos):
    """
    :param tree_path: List of indices representing the path in the instrument tree.
    :param instrument_infos: String containing information about instrument.
    :return: Updated tree path.
    This method searches for an ionisation by iterating through the keys in the instrument tree specified by the tree_path.
    If a key is found in the instrument_infos string, it is appended to the tree_path and returned.
    If no matching key is found, 'unknown' is appended to the tree_path and returned.
    """
    try:
        # Iterating over keys in the instrument tree at the specified path
        for key in globals_vars.instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]].keys():
            # If the key is present in instrument_infos, append it to the path and return
            if re.search(rf"(\b|^|$){key}(\b|^|$)", instrument_infos):
                tree_path.append(key)
                return tree_path
        # If no match found, append 'unknown' to the path and return
        tree_path.append('unknown')
        return tree_path
    except:
        # If anything goes wrong, return None
        return None

def make_tree_path(instrument_infos):
    """
    :param instrument_infos: A list of dictionaries containing information about instruments
    :return: A list representing the tree path based on the given instrument information

    This method takes in a list of dictionaries containing information about instruments. It returns a list representing the tree path based on the given instrument information. The tree
    * path is constructed by searching for specific attributes within the instrument information and appending them to the tree path list.

    The method first initializes an empty list called tree_path. Then, it calls various search functions to populate the tree_path in the following order:

    1. search_for_brand: searches for the brand attribute in the instrument information and appends it to the tree_path
    2. search_for_model: searches for the model attribute in the instrument information and appends it to the tree_path
    3. search_for_spectrum_type: searches for the spectrum type attribute in the instrument information and appends it to the tree_path
    4. search_for_instrument_type: searches for the instrument type attribute in the instrument information and appends it to the tree_path
    5. search_for_ionisation: searches for the ionisation attribute in the instrument information and appends it to the tree_path

    Finally, the method returns the tree_path list.

    Example usage:

    instrument_info = [
        {"brand": "Brand1", "model": "Model1", "spectrum_type": "Type1", "instrument_type": "TypeA", "ionisation": "IonA"},
        {"brand": "Brand2", "model": "Model2", "spectrum_type": "Type2", "instrument_type": "TypeB", "ionisation": "IonB"},
    ]

    tree = make_tree_path(instrument_info)
    print(tree)
    # Output: ["Brand1", "Model1", "Type1", "TypeA", "IonA", "Brand2", "Model2", "Type2", "TypeB", "IonB"]
    """
    # Initialize an empty list for the tree path
    tree_path = []

    # Search for the brand attribute and append it to the tree path
    tree_path = search_for_brand(tree_path, instrument_infos)

    # If not found, return None
    if not tree_path:
        return None

    # Search for the model attribute and append it to the tree path
    tree_path = search_for_model(tree_path, instrument_infos)

    # If not found, return None
    if not tree_path:
        return None

    # Search for the spectrum type attribute and append it to the tree path
    tree_path = search_for_spectrum_type(tree_path, instrument_infos)

    # If not found, return None
    if not tree_path:
        return None

    # Search for the instrument type attribute and append it to the tree path
    tree_path = search_for_instrument_type(tree_path, instrument_infos)

    # If not found, return None
    if not tree_path:
        return None

    # Search for the ionisation attribute and append it to the tree path
    tree_path = search_for_ionisation(tree_path, instrument_infos)

    # If not found, return the tree path as it is
    if not tree_path:
        return None

    # Return the constructed tree path
    return tree_path

def normalize_instruments_and_resolution(metadata_dict):
    """
    This function retrieves the instrument and resolution information from the provided
    metadata dictionary, modifies them and updates the dictionary with the adjusted values.

    :param metadata_dict: A dictionary containing metadata information. It should have
    the keys "INSTRUMENT", "INSTRUMENTTYPE", and "RESOLUTION". These are the keys which will be
    modified by this method.
    :return: The modified metadata_dict dictionary.
    """

    # Get the cleaned instrument information for the given metadata
    instrument_infos = clean_spectrum_instrument_info(metadata_dict)

    # Ensur academia usage of the instrument information (additional information between '.')
    instrument_infos = f". {instrument_infos} ."

    # Generate path in the instrument catalogue by using instrument information
    tree_path = make_tree_path(instrument_infos)

    # If no path was found in the catalogue, return the metadata unaltered
    if not tree_path:
        return metadata_dict

    # Try to retrieve the instrument's resolution and specific solution
    try:
        # Determine if the instrument has high or low resolution, if neither can be determined, mark as unknown
        resolution = "high" if "high" in globals_vars.instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]][tree_path[4]] else "low" if "low" in globals_vars.instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]][tree_path[4]] else "unknown"

        # Retrieve the corresponding solution for the determined resolution
        solution = globals_vars.instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]][tree_path[4]][resolution]["SOLUTION"]

        # Split the solution string into separate components
        solution = solution.split(',')
    except:
        # If an error occurred during retrieval, return the metadata unaltered
        return metadata_dict

    # Update the metadata dictionary with the retrieved solution
    metadata_dict["INSTRUMENT"] = solution[0].strip()
    metadata_dict["INSTRUMENTTYPE"] = solution[1].strip()
    metadata_dict["RESOLUTION"] = solution[2].strip()

    # If the instrument type contains the ionisation method, separate it and place it in its own field
    if len(solution[1].split('-')) >= 2:
        metadata_dict["IONIZATION"] = solution[1].split('-')[1].strip()

    # Return the updated metadata dictionary
    return metadata_dict
