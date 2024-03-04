from .keys_convertor import *
import concurrent.futures
from tqdm import tqdm
import re

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*:.[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:\n]*?):\s*([^:\n]*)(?:\n|$)")

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
    metadata,peak_list = None,None

    match = re.search(metadata_peak_list_split_pattern, spectrum)

    if match:
        metadata, peak_list = match.group(1), match.group(2)

    return metadata,peak_list

def check_for_metadata_in_comments(metadata_matches):
    """
    :param metadata_matches: A list of metadata matches to be checked.
    :return: A modified list of metadata matches if any matches were found in the comments, otherwise False.

    This method checks if the given metadata matches contain metadata found in comments. It iterates over each match and checks if a comment field exists. If a comment field exists, the
    * match is added to a new list, and if the field contains a "=" character, sub-fields are checked for matches as well. If sub-fields are found, they are added to the new list as separate
    * tuples. If no metadata is found in the comments, False is returned.

    Example:
    metadata_matches = [('Comment', 'Field=Some value'), ('Another Comment', 'Field=Value')]
    check_for_metadata_in_comments(metadata_matches)

    Output:
    [('Comment', 'Field=Some value'), ('Another Comment', 'Field=Value')]

    """
    new_metadata_matches = []
    # Check if comment filed exist
    for match in metadata_matches:
        if re.search(comment_pattern, match[0]):
            new_metadata_matches.append(match)
            if "=" in match[1]:
                sub_fields_matches = re.findall(sub_fields_pattern, match[1])
                if sub_fields_matches:
                    for sub_fields_match in sub_fields_matches:
                        non_empty_tuple = tuple(group for group in sub_fields_match if group)
                        if len(non_empty_tuple) == 2:
                            new_metadata_matches.append(non_empty_tuple)
                else:
                    return False
            else:
                new_metadata_matches.append(match)
        else:
            new_metadata_matches.append(match)

    return new_metadata_matches if new_metadata_matches else False

def metadata_to_dict(metadata):
    """
    Converts metadata string to a dictionary.

    :param metadata: The metadata string to be converted.
    :return: A dictionary representation of the metadata.

    """
    metadata_dict = {}

    metadata_matches = re.findall(metadata_pattern,metadata)

    if metadata_matches:
        temp = check_for_metadata_in_comments(metadata_matches)
        temp = [t for t in temp if not re.search(computed_pattern, t[0])]
        if temp != False:
            metadata_matches = temp

        for match in metadata_matches:
            metadata_dict[re.sub(metadata_fields_name_pattern, '', match[0]).lower().strip()] = re.sub(metadata_strip_value_pattern,"",match[1])

        return metadata_dict

    return metadata_dict

def peak_list_to_array(peak_list):
    """
    :param peak_list: A string representing a list of peaks. Each peak in the list is represented by a pair of numbers (x, y), where x and y are floats.
    :return: A list of peak arrays. Each peak array is a pair of floats (x, y).

    """
    peaks_match = re.findall(peak_list_split_pattern, peak_list)

    if peaks_match:
        peak_array = [[float(i), float(j)] for i, j in peaks_match]
        return peak_array
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
    if not metadata or not peak_list:
        return {},[]
    else:
        metadata_dict = metadata_to_dict(metadata)
        if not metadata_dict:
            return {},[]
        else:
            peak_list = peak_list_to_array(peak_list)
            if peak_list:
                return metadata_dict, peak_list
            else:
                return {},[]

def msp_to_json(spectrum):
    """
    Convert an MSP spectrum to JSON format.

    :param spectrum: The input MSP spectrum.
    :return: The JSON representation of the spectrum.
    """
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    metadata, peak_list = structure_metadata_and_peak_list(metadata, peak_list)
    if not metadata or not peak_list:
        return None
    metadata["peaks"] = peak_list

    # metadata = convert_keys(metadata)

    return metadata.keys()

def msp_to_json_processing(FINAL_MSP):
    """
    Process a list of spectrums and convert them to JSON format.

    :param spectrum_list: List of spectrums to be processed.
    :return: List of spectrums in JSON format.
    """

    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(FINAL_MSP), unit=" spectrums", colour="green", desc="{:>70}".format("converting MSP spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(FINAL_MSP), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_MSP[i:i + chunk_size]
            results = list(executor.map(msp_to_json, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final