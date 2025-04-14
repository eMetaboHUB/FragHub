from scripts.calculate_maximized_chunk_size import *
from scripts.convertors.keys_convertor import *
import concurrent.futures
import scripts.globals_vars
import re


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
    match = re.search(scripts.globals_vars.metadata_peak_list_split_pattern_mgf, spectrum)
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
    metadata_matches = re.findall(scripts.globals_vars.metadata_pattern_mgf, metadata)
    del metadata

    # Check if there are matches
    if metadata_matches:
        # Iterate over each match
        for match in metadata_matches:
            # Remove the field name from the match, convert to lower case and strip spaces, then assign to dictionary
            # The assumed structure of each match is (field name, value) pair
            metadata_dict[re.sub(scripts.globals_vars.metadata_fields_name_pattern, '', match[0]).lower().strip()] = re.sub(scripts.globals_vars.metadata_strip_value_pattern, "", match[1])

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
    peaks_match = re.findall(scripts.globals_vars.peak_list_split_pattern, peak_list)
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

def mgf_to_dict(spectrum):
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


def mgf_to_dict_processing(FINAL_MGF, progress_callback=None, total_items_callback=None, prefix_callback=None,
                           item_type_callback=None):
    """
    Converts MGF spectrums to JSON format, with support for progress reporting via callbacks.
    :param FINAL_MGF: A list of MGF spectrums to be converted.
    :param progress_callback: A function to update the progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: A list of JSON-formatted spectrums.
    """
    start = 0
    end = len(FINAL_MGF)

    # Set total items via callback if provided
    if total_items_callback:
        total_items_callback(end, 0)  # total = end, completed = 0

    # Update the prefix dynamically via callback if provided
    if prefix_callback:
        prefix_callback("Parsing MGF spectrums:")

    # Specify the type of items being processed via callback if provided
    if item_type_callback:
        item_type_callback("spectra")

    # Variable to track progress
    processed_items = 0

    # Calculate the chunk size once, outside the loop
    chunk_size = calculate_maximized_chunk_size(data_list=FINAL_MGF)

    # Process the data in chunks
    while start < end:
        # Use ThreadPoolExecutor to process the chunk
        with concurrent.futures.ThreadPoolExecutor() as executor:
            FINAL_MGF[start:start + chunk_size] = list(executor.map(mgf_to_dict, FINAL_MGF[start:start + chunk_size]))

        # Filter out None results
        FINAL_MGF[start:start + chunk_size] = [item for item in FINAL_MGF[start:start + chunk_size] if item is not None]

        # Update progress via callback if provided
        processed_items += min(chunk_size, end - start)  # Ensure not to exceed the total size
        if progress_callback:
            progress_callback(processed_items)

        # Move to the next chunk
        start += chunk_size

    return FINAL_MGF


