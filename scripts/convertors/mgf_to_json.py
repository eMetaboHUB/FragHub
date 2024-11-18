from .keys_convertor import *
import concurrent.futures
from tqdm import tqdm
import re

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*=.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:\n]*?)=\s*([^\n]*)(?:\n|$)")

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def extract_metadata_and_peak_list(spectrum):
    """
    The function inputs a spectrum and extracts the metadata and peak list out of it

    :param spectrum: The input is a spectrum from which the function extracts metadata and peak list

    :return: The function returns a tuple containing the metadata and peak list. If no match is found, the function
             returns None for both
    """
    # The initial metadata and peak list are set as None
    metadata, peak_list = None, None

    # The function then searches for a match in the spectrum based on the regex metadata_peak_list_split_pattern
    match = re.search(metadata_peak_list_split_pattern, spectrum)
    del spectrum

    # If the search finds a match, metadata and peak list are defined as group 1 & 2 of the match respectively
    if match:
        metadata, peak_list = match.group(1), match.group(2)

    # The results are returned as a tuple (metadata, peak_list)
    return metadata, peak_list

def metadata_to_dict(metadata):
    """
    Convert metadata to a dictionary.
    :param metadata: The metadata string.
    :type metadata: str
    :return: The metadata as a dictionary.
    :rtype: dict
    """
    # Initialize an empty dictionary to hold the metadata
    metadata_dict = {}

    # Use regular expression to find all metadata matches in the given string
    metadata_matches = re.findall(metadata_pattern, metadata)
    del metadata

    # Check if there are matches
    if metadata_matches:
        # Iterate over each match
        for match in metadata_matches:
            # Remove the field name from the match, convert to lower case and strip spaces, then assign to dictionary
            # The assumed structure of each match is (field name, value) pair
            metadata_dict[re.sub(metadata_fields_name_pattern, '', match[0]).lower().strip()] = re.sub(metadata_strip_value_pattern, "", match[1])

        # Return the dictionary with metadata
        return metadata_dict

    # In case there are no matches, return an empty dictionary
    return metadata_dict

def peak_list_to_array(peak_list):
    """
    Convert a peak list to an array.
    :param peak_list: A string containing a list of peaks in the format "x1,y1 x2,y2 x3,y3 ...".
    :return: A list of peak coordinates in the format [[x1, y1], [x2, y2], [x3, y3], ...]. If the peak list is empty or
             doesn't match the expected format, an empty list is returned.
    """

    # Use regular expression to find all peak pairs in the given string
    peaks_match = re.findall(peak_list_split_pattern, peak_list)
    del peak_list

    # If matches exist
    if peaks_match:
        # Convert the matches into an array of lists with x, y as float
        peak_array = [[float(i), float(j)] for i, j in peaks_match]
        del peaks_match
        return peak_array
    else:
        # If no matches found or the peak list is empty, return an empty list
        return []

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    :param metadata: The metadata for the structure.
    :param peak_list: The list of peaks for the structure.
    :return: A tuple containing a dictionary representing the structured metadata and a list representing the structured peak list. If either the metadata or the peak list is empty, an empty
    * dictionary and an empty list will be returned, respectively.
    """

    # If metadata or peak list is empty return empty dictionary and list respectively
    if not metadata or not peak_list:
        return {}, []
    else:
        # Convert the metadata into a dictionary
        metadata_dict = metadata_to_dict(metadata)
        del metadata
        # If metadata dictionary is empty return empty dictionary and list
        if not metadata_dict:
            return {}, []
        else:
            # Convert peak list to array
            peak_list = peak_list_to_array(peak_list)
            # If peak list exist return both metadata dictionary and peak list
            if peak_list:
                return metadata_dict, peak_list
            else:
                # If peak list does not exist return empty dictionary and list
                return {}, []

def mgf_to_json(spectrum):
    """
    Convert MGF spectrum to JSON format.
    :param spectrum: The MGF spectrum to convert.
    :type spectrum: dict
    :return: The converted spectrum in JSON format, or None if conversion fails.
    :rtype: dict or None
    """

    # Extract the metadata and peak list from the spectrum
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    del spectrum

    # Structure the extracted metadata and peak list
    metadata, peak_list = structure_metadata_and_peak_list(metadata, peak_list)

    # If the metadata or peak list is None or not present, return None
    if not metadata or not peak_list:
        return None

    # Add the peak list to the metadata dictionary under the key 'peaks'
    metadata["peaks"] = peak_list

    # Convert the keys in the metadata to a standardized format
    metadata = convert_keys(metadata)

    # Return the modified and structured metadata
    return metadata

def mgf_to_json_processing(FINAL_MGF):
    """
    Converts MGF spectrums to JSON format.
    :param FINAL_MGF: A list of MGF spectrums to be converted.
    :return: A list of JSON-formatted spectrums.
    """
    start = 0
    end = len(FINAL_MGF)

    # Initialize a progress bar to keep track of and visualize the conversion process
    progress_bar = tqdm(total=end, unit=" spectrums", colour="green", desc="{:>70}".format("converting MGF spectrums"))

    while start < end:
        chunk_size = min(end - start, 5000)

        # Use ThreadPoolExecutor to process the chunk
        with concurrent.futures.ThreadPoolExecutor() as executor:
            FINAL_MGF[start:start + chunk_size] = list(executor.map(mgf_to_json, FINAL_MGF[start:start + chunk_size]))

        # Filter out None results
        FINAL_MGF[start:start + chunk_size] = [item for item in FINAL_MGF[start:start + chunk_size] if item is not None]

        # Update progress bar
        progress_bar.update(chunk_size)

        # Move to the next chunk
        start += chunk_size

    progress_bar.close()
    return FINAL_MGF
