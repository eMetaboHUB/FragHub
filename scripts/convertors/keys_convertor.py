import pandas as pd
import globals_vars
import os


def convert_keys(metadata_dict):
    """
    Convert keys in metadata_dict based on the provided keys_dict and keys_list.
    Dictionaries are a built-in data type in Python used to store collections of items that are mapped to keys.

    :param metadata_dict: A dictionary containing metadata information.

    :return: A dictionary with converted keys based on the predefined keys_dict and keys_list.
    """
    # Creating a dictionary comprehension that converts all keys from metadata to lower case
    # and matches them with the keys available in keys_dict and keys_list.
    # The output is a dictionary where the keys are mapped from keys_dict and the values come from metadata_dict.
    converted = {globals_vars.keys_dict[key.lower()]: val for key, val in metadata_dict.items() if key.lower() in keys_dict and keys_dict[key.lower()] in keys_list}
    del metadata_dict

    # After initial conversion, there might still be some keys from keys_list that are not in the converted dictionary.
    # This line adds those missing keys to the converted dictionary with an empty string ("") as their value.
    converted.update({key: "" for key in globals_vars.keys_list if key not in converted})

    # Returning the result that is a dictionary with matching and new keys.
    return converted
