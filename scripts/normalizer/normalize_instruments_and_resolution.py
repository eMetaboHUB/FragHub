import json
import re

global instrument_tree
with open('../datas/instruments_tree.json', 'r') as f:
    instrument_tree = json.load(f)

def clean_instrument(instrument):
    """
    Clean the instrument name string by removing specific prefixes and modifying the format.

    :param instrument: The instrument name string to be cleaned.
    :return: The cleaned instrument name string.
    """
    instrument = re.sub("-tof", "tof", instrument)
    instrument = re.sub("q-", "q", instrument)
    instrument = re.sub("q exactive", "qexactive", instrument)
    instrument = re.sub("applied biosystems", "ab sciex", instrument)
    instrument = re.sub(" ab ", " ab sciex ", instrument)
    instrument = re.sub("sciex", " ab sciex ", instrument)
    instrument = re.sub("triple(-| )?tof", "qqq", instrument)
    instrument = re.sub("triple(-| )?quad", "qqq", instrument)
    instrument = re.sub("... uplc ...", "", instrument)

    return instrument

def clean_instrument_type(instrument_type):
    """
    Cleans the given instrument string by removing specific substrings and replacing hyphens with spaces.

    :param instrument_str: The instrument string to be cleaned.
    :type instrument_str: str
    :return: The cleaned instrument string.
    :rtype: str
    """
    instrument_type = re.sub("-tof", "tof", instrument_type)
    instrument_type = re.sub("q-", "q", instrument_type)
    instrument_type = re.sub("-", " ", instrument_type)
    instrument_type = re.sub("q exactive", "qexactive", instrument_type)
    instrument_type = re.sub("applied biosystems", "ab sciex", instrument_type)
    instrument_type = re.sub(" ab ", " ab sciex ", instrument_type)
    instrument_type = re.sub("sciex", " ab sciex ", instrument_type)
    instrument_type = re.sub("triple(-| )?tof", "qqq", instrument_type)
    instrument_type = re.sub("triple(-| )?quad", "qqq", instrument_type)
    instrument_type = re.sub("... uplc ...", "", instrument_type)

    return instrument_type

def clean_spectrum_instrument_info(metadata_dict):
    """
    Cleans the spectrum instrument information from the given metadata dictionary.

    :param metadata_dict: The dictionary containing the metadata information.
    :return: The cleaned instrument information.
    """
    instrument = metadata_dict['INSTRUMENT'].lower()
    instrument_type = metadata_dict["INSTRUMENTTYPE"].lower()

    instrument = clean_instrument(instrument)
    instrument_type = clean_instrument_type(instrument_type)


    instrument_infos = instrument + " " + instrument_type
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
    for key in instrument_tree.keys():
        if f" {key} " in instrument_infos:
            tree_path.append(key)
            return tree_path

    tree_path.append('unknown')
    return tree_path

def search_for_model(tree_path, instrument_infos):
    """
    Search for a model in the instrument tree based on the given tree path and instrument infos.

    :param tree_path: The path to the model in the instrument tree.
    :param instrument_infos: The instrument infos used to search for the model.
    :return: The updated tree path with the found model or 'unknown' if not found.
    """
    for key in instrument_tree[tree_path[0]].keys():
        if f" {key} " in instrument_infos:
            tree_path.append(key)
            return tree_path

    tree_path.append('unknown')
    return tree_path

def search_for_spectrum_type(tree_path, instrument_infos):
    """
    Method to search for spectrum type in a given tree path and instrument infos.

    :param tree_path: List of keys representing the tree path.
    :param instrument_infos: String containing instrument information.
    :return: Updated tree path with found spectrum type, or 'unknown' if not found.
    """
    for key in instrument_tree[tree_path[0]][tree_path[1]].keys():
        if f" {key} " in instrument_infos:
            tree_path.append(key)
            return tree_path

    tree_path.append('unknown')
    return tree_path

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
    for key in instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]].keys():
        if f" {key} " in instrument_infos:
            tree_path.append(key)
            return tree_path

    tree_path.append('unknown')
    return tree_path

def search_for_ionisation(tree_path, instrument_infos):
    """
    :param tree_path: List of indices representing the path in the instrument tree.
    :param instrument_infos: String containing information about instrument.
    :return: Updated tree path.

    This method searches for an ionisation by iterating through the keys in the instrument tree specified by the tree_path.
    If a key is found in the instrument_infos string, it is appended to the tree_path and returned.
    If no matching key is found, 'unknown' is appended to the tree_path and returned.
    """
    for key in instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]].keys():
        if f" {key} " in instrument_infos:
            tree_path.append(key)
            return tree_path

    tree_path.append('unknown')
    return tree_path

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
    tree_path = []

    tree_path = search_for_brand(tree_path, instrument_infos)
    tree_path = search_for_model(tree_path, instrument_infos)
    tree_path = search_for_spectrum_type(tree_path, instrument_infos)
    tree_path = search_for_instrument_type(tree_path, instrument_infos)
    tree_path = search_for_ionisation(tree_path, instrument_infos)

    return tree_path

def normalize_instruments_and_resolution(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
        - The dictionary should have the keys "INSTRUMENT", "INSTRUMENTTYPE", and "RESOLUTION".
        - The values of these keys will be modified by this method.
    :return: The modified metadata_dict dictionary.

    This method normalizes the instrument and resolution information in the given metadata_dict. It retrieves the instrument information from the metadata_dict, cleans it, and updates the
    * relevant keys in the metadata_dict with the cleaned values.
    """
    instrument_infos = clean_spectrum_instrument_info(metadata_dict)

    tree_path = make_tree_path(instrument_infos)

    solution = instrument_tree[tree_path[0]][tree_path[1]][tree_path[2]][tree_path[3]][tree_path[4]]["SOLUTION"]
    solution = solution.split(',')

    metadata_dict["INSTRUMENT"] = solution[0].strip()
    metadata_dict["INSTRUMENTTYPE"] = solution[1].strip()
    metadata_dict["RESOLUTION"] = solution[2].strip()

    return metadata_dict