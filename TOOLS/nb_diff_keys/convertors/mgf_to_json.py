from .keys_convertor import *
import concurrent.futures
from tqdm import tqdm
import re

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*=.[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:\n]*?)=\s*([^:\n]*)(?:\n|$)")

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def extract_metadata_and_peak_list(spectrum):
    """
    Extracts metadata and peak list from the given spectrum.

    :param spectrum: The input spectrum.
    :return: A tuple containing metadata and peak list. Returns None for both if no match is found.
    """
    metadata, peak_list = None, None

    match = re.search(metadata_peak_list_split_pattern, spectrum)

    if match:
        metadata, peak_list = match.group(1), match.group(2)

    return metadata, peak_list

def metadata_to_dict(metadata):
    """
    Convert metadata to a dictionary.

    :param metadata: The metadata string.
    :type metadata: str
    :return: The metadata as a dictionary.
    :rtype: dict
    """
    metadata_dict = {}

    metadata_matches = re.findall(metadata_pattern, metadata)

    if metadata_matches:
        for match in metadata_matches:
            metadata_dict[re.sub(metadata_fields_name_pattern, '', match[0]).lower().strip()] = re.sub(metadata_strip_value_pattern, "", match[1])

        return metadata_dict

    return metadata_dict

def peak_list_to_array(peak_list):
    """
    Convert a peak list to an array.

    :param peak_list: A string containing a list of peaks in the format "x1,y1 x2,y2 x3,y3 ...".
    :return: A list of peak coordinates in the format [[x1, y1], [x2, y2], [x3, y3], ...]. If the peak list is empty or
             doesn't match the expected format, an empty list is returned.
    """
    peaks_match = re.findall(peak_list_split_pattern, peak_list)

    if peaks_match:
        peak_array = [[float(i), float(j)] for i, j in peaks_match]
        return peak_array
    else:
        return []

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    :param metadata: The metadata for the structure.
    :param peak_list: The list of peaks for the structure.
    :return: A tuple containing a dictionary representing the structured metadata and a list representing the structured peak list. If either the metadata or the peak list is empty, an empty
    * dictionary and an empty list will be returned, respectively.
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

def mgf_to_json(spectrum):
    """
    Convert MGF spectrum to JSON format.

    :param spectrum: The MGF spectrum to convert.
    :type spectrum: dict
    :return: The converted spectrum in JSON format, or None if conversion fails.
    :rtype: dict or None
    """
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    metadata, peak_list = structure_metadata_and_peak_list(metadata, peak_list)
    if not metadata or not peak_list:
        return None
    metadata["peaks"] = peak_list
    # metadata = convert_keys(metadata)

    return metadata.keys()

def mgf_to_json_processing(FINAL_MGF):
    """
    Converts MGF spectrums to JSON format.

    :param FINAL_MGF: A list of MGF spectrums to be converted.
    :return: A list of JSON-formatted spectrums.

    Example usage:
    ```
    final_json = mgf_to_json_processing(FINAL_MGF)
    ```
    """
    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(FINAL_MGF), unit=" spectrums", colour="green", desc="{:>70}".format("converting MGF spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(FINAL_MGF), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_MGF[i:i + chunk_size]
            results = list(executor.map(mgf_to_json, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final