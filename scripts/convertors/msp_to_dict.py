from calculate_maximized_chunk_size import *
from .keys_convertor import *
import concurrent.futures
import re

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*:.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:]*):(?: )?([^\n]*)(?:\n|$)")

global computed_pattern
computed_pattern = re.compile(r"computed", flags=re.IGNORECASE)

global comment_pattern
comment_pattern = re.compile(r'comment.*', flags=re.IGNORECASE)

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global sub_fields_pattern
sub_fields_pattern = re.compile(r"(\S+?)=\"([^\"]*)\"|\"(\w+?)=([^\"]*)\"|\"([^\"]*?)=([^\"]*)\"|(\S+?)=(\d+(?:[.,]\d*)?)|(\S+?)=(.*?)(?:;|\n|$)")

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

def extract_metadata_and_peak_list(spectrum):
    """
    Extracts the metadata and peak list from the given spectrum.
    :param spectrum: The spectrum from which to extract the metadata and peak list.
    :type spectrum: str
    :return: The metadata and peak list extracted from the spectrum. If the regex pattern does not match, both metadata and peak list will be None.
    :rtype: tuple(str, str)
    """

    # Initialize the metadata and peak_list as None
    metadata, peak_list = None, None

    # Use a regex pattern to search for matches within the input spectrum
    match = re.search(metadata_peak_list_split_pattern, spectrum)

    # If a match is found, assign the groups to metadata and peak_list variables
    if match:
        metadata, peak_list = match.group(1), match.group(2)

    # Return the metadata and peak_list
    return metadata, peak_list

def check_for_metadata_in_comments(metadata_matches):
    """
    This function checks for metadata within comments fields and collects them.

    :param metadata_matches: A list of metadata matches that needs to be checked.
    :return: A modified list of metadata matches if any matches were found within the comments.
             If no matches are found, the function returns False.
    """

    new_metadata_matches = []  # List to store modified/filtered metadata matches

    # Iterating over each metadata match
    for match in metadata_matches:
        # If a comment field exists in the given match
        if re.search(comment_pattern, match[0]):
            # Add the match to the new list
            new_metadata_matches.append(match)

            # Check if the comment field contains any "=" character
            if "=" in match[1]:
                # Look for sub-fields in the comment
                sub_fields_matches = re.findall(sub_fields_pattern, match[1])

                # If sub-fields are found, add them in the new list as separate tuples
                if sub_fields_matches:
                    for sub_fields_match in sub_fields_matches:
                        # Filter to keep non-empty groups only
                        non_empty_tuple = tuple(group for group in sub_fields_match if group)
                        if len(non_empty_tuple) == 2:
                            new_metadata_matches.append(non_empty_tuple)
                else:
                    # If no subfield matches are found, return False
                    return False
            else:
                # If no "=" character is found, simply add the match to the list
                new_metadata_matches.append(match)
        else:
            # If comment field doesn't exist, add the match to the list
            new_metadata_matches.append(match)

    # Return the modified metadata matches if any found, otherwise return False
    return new_metadata_matches if new_metadata_matches else False

def metadata_to_dict(metadata):
    """
    Converts metadata string to a dictionary. This function takes a metadata string
    as input and returns a dictionary. The metadata string is expected to follow a specific
    pattern which this function will try to match and use to generate the dictionary

    :param metadata: The string representing the metadata
    :return: A dictionary where keys represent metadata field names and values are the respective data
    """

    # Initializing an empty dictionary to store our metadata in key-value pairs
    metadata_dict = {}

    # Using regex to find all matches in the metadata string that follow the specified pattern
    metadata_matches = re.findall(metadata_pattern, metadata)

    # If we have found matches, further process them. If not, return an empty dictionary
    if metadata_matches:

        # Checking for metadata in comments and returning a new list with only valid metadata
        temp = check_for_metadata_in_comments(metadata_matches)

        # Removing any computational metadata as they do not provide useful information
        if temp:
            temp = [t for t in temp if not re.search(computed_pattern, t[0])]

        # If temp isn't False then it contains parsed metadata for processing
        if temp != False:
            # Reassign metadata_matches to hold only valid matches for further processing
            metadata_matches = temp

        # Processing each match found and preparing a dictionary
        for match in metadata_matches:
            # Extracting field name from the match using sub-regex, also making it lowercase and removing starting and trailing spaces
            field_name = re.sub(metadata_fields_name_pattern, '', match[0]).lower().strip()

            # Extracting field value from the match using sub-regex
            field_value = re.sub(metadata_strip_value_pattern, "", match[1])

            # Storing the field name and value in our dictionary
            metadata_dict[field_name] = field_value

        # Returning the dictionary
        return metadata_dict

    # If there are no matches found, return an empty dictionary
    return metadata_dict

def peak_list_to_array(peak_list):
    """
    :param peak_list: A string representing a list of peaks. Each peak in the list is represented by a pair of numbers (x, y), where x and y are floats.
    :return: A list of peak arrays. Each peak array is a pair of floats (x, y).

    The function takes a 'peak_list' argument which represents a list of peaks stored as a string. The objective is to convert this string representation of peak list into an array.

    The function uses regular expressions to find all peaks in the string that match the defined pattern. If the function finds matches, the list of peaks is converted into an array.

    Each peak in the list initially represented in the form of string (x, y), is converted into a list of two floats [x, y].

    The resulting list of peaks is then returned.

    If there are no matches found (meaning no peaks were found in the peak_list string), the function returns an empty list.
    """
    peak_list = re.findall(peak_list_split_pattern, peak_list)  # use regex to find all peaks in the string

    # If peak matches found return list of lists, where inner lists are pair of floats
    if peak_list:
        peak_array = [[float(i), float(j)] for i, j in peak_list]
        return peak_array
    # If peak matches not found return an empty list
    else:
        return []

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    :param metadata: The metadata to be structured and converted to a dictionary.
    :param peak_list: The peak list to be structured and converted to an array.
    :return: A tuple containing the structured metadata as a dictionary and the structured peak list as an array.

    The function `structure_metadata_and_peak_list` takes in two parameters `metadata` and `peak_list`, and returns a tuple containing the structured metadata and peak list.

    To structure the metadata, the function first checks if either the `metadata` parameter or the `peak_list` parameter is empty. If any of them is empty, the function returns an empty
    * dictionary and an empty array.

    Then, the function converts the `metadata` parameter to a dictionary using the `metadata_to_dict` function.

    If the `metadata_dict` is not empty, the function proceeds to convert the `peak_list` parameter to an array using the `peak_list_to_array` function.

    If the conversion is successful and the `peak_list` is not empty, the function returns the structured metadata dictionary and the structured peak list array.

    If the conversion is unsuccessful or the `peak_list` is empty, the function returns an empty dictionary and an empty array.
    """
    # Check if the metadata or peak list is empty
    if not metadata or not peak_list:
        # If either is empty, return empty dictionary and list
        return {}, []
    else:
        # Structure the metadata into dictionary form
        metadata = metadata_to_dict(metadata)
        # Check if the new structured metadata is empty
        if not metadata:
            # If it's empty, return an empty dictionary and list
            return {}, []
        else:
            # If it's not empty, structure the peak list into an array.
            peak_list = peak_list_to_array(peak_list)
            # If the structured peak list is not empty, return both structured metadata and peak list.
            if peak_list:
                return metadata, peak_list
            else:
                # If the peak list is, nonetheless, empty, return an empty dictionary and list.
                return {}, []

def msp_to_dict(spectrum):
    # Function to convert a given MSP spectrum to its corresponding JSON format

    """
    This function takes in a MSP spectrum, extracts its metadata and peak list, structures them,
    and if valid metadata and peak list exist, saves the peak list into the metadata and converts
    the keys in metadata to canonical forms, and finally returns the processed metadata,
    which is a JSON-like structured Python dictionary.

    :param spectrum: The input MSP spectrum.
    :return: The JSON-like structured Python dictionary representation of the spectrum, or None if no valid metadata or peak list exists.
    """

    # Extract general (meta) info and raw peak info from the MSP spectrum
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    del spectrum

    # Convert the unstructured metadata and peak list extracted from the spectrum to structured forms
    metadata, peak_list = structure_metadata_and_peak_list(metadata, peak_list)

    # Check if valid metadata and peak list exist
    # If not, return None; if yes, continue to the next steps
    if not metadata or not peak_list:
        return None

    # Save the peak list into the metadata dictionary
    metadata["peaks"] = peak_list

    # Convert the keys in the metadata dictionary to canonical (standard) forms
    metadata = convert_keys(metadata)

    # Return the processed metadata, which is a JSON-like structured Python dictionary
    return metadata

def msp_to_dict_processing(FINAL_MSP, progress_callback=None, total_items_callback=None, prefix_callback=None,
                           item_type_callback=None):
    """
    Process a list of MSP spectrums and convert them to JSON format, with support for progress callbacks.
    :param FINAL_MSP: List of MSP spectrums to be processed.
    :param progress_callback: A function to update the progress.
    :param total_items_callback: A function to set the total number of items.
    :param prefix_callback: A function to dynamically set the prefix for the operation.
    :param item_type_callback: A function to specify the type of items processed (e.g., "spectra").
    :return: List of spectrums in JSON format.
    """
    start = 0
    end = len(FINAL_MSP)

    # Initialiser le total des items via callback (si défini)
    if total_items_callback:
        total_items_callback(end, 0)  # total = end, completed = 0

    # Mise à jour dynamique du préfixe (si défini)
    if prefix_callback:
        prefix_callback("Parsing MSP spectrums:")

    # Mise à jour du type d'éléments (si défini)
    if item_type_callback:
        item_type_callback("spectra")

    # Initialisation de la variable pour suivre la progression
    processed_items = 0

    # Calculer la taille des chunks une seule fois au départ
    chunk_size = calculate_maximized_chunk_size(data_list=FINAL_MSP)

    # Traitement en morceaux (chunks) pour les spectres MSP
    while start < end:
        # Utilisation de ThreadPoolExecutor pour le traitement parallèle
        with concurrent.futures.ThreadPoolExecutor() as executor:
            FINAL_MSP[start:start + chunk_size] = list(executor.map(msp_to_dict, FINAL_MSP[start:start + chunk_size]))

        # Filtrer les résultats None
        FINAL_MSP[start:start + chunk_size] = [item for item in FINAL_MSP[start:start + chunk_size] if item is not None]

        # Mise à jour de la progression
        processed_items += min(chunk_size, end - start)  # S'assurer de ne pas dépasser la taille totale
        if progress_callback:
            progress_callback(processed_items)

        # Passer au prochain morceau
        start += chunk_size

    return FINAL_MSP


