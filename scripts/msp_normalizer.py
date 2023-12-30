from values_normalizer import *
import concurrent.futures
from Filters_NEW import *
from tqdm import tqdm
import pandas as pd
import numpy as np
import time
import os
import re

np.set_printoptions(suppress=True)

global keys_dict
Key_dataframe = pd.read_csv(os.path.abspath("../datas/key_to_convert.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
keys_dict = dict(zip(Key_dataframe['known_synonym'], Key_dataframe['fraghub_default'].str.upper()))

global keys_list
keys_list = ['FILENAME',
             'PREDICTED',
             'FRAGHUBID',
             'SPECTRUMID',
             'RESOLUTION',
             'SYNON',
             'CHARGE',
             'IONIZATION',
             'MSLEVEL',
             'FRAGMENTATIONMODE',
             'NAME',
             'PRECURSORMZ',
             'EXACTMASS',
             'AVERAGEMASS',
             'PRECURSORTYPE',
             'INSTRUMENTTYPE',
             'INSTRUMENT',
             'SMILES',
             'INCHI',
             'INCHIKEY',
             'COLLISIONENERGY',
             'FORMULA',
             'RETENTIONTIME',
             'IONMODE',
             'COMMENT',
             'NUM PEAKS']

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile("([\s\S]*:.[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global sub_fields_pattern
sub_fields_pattern = re.compile("(\S+?)=\"([^\"]*)\"|\"(\w+?)=([^\"]*)\"|\"([^\"]*?)=([^\"]*)\"|(\S+?)=(\d+(?:[.,]\d*)?)|(\S+?)=(.*?)(?:;|\n|$)")

global metadata_pattern
metadata_pattern = re.compile("([^:\n]*?):\s*([^:\n]*)(?:\n|$)")

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile("^\"|\"$")

global peak_list_split_pattern
peak_list_split_pattern = re.compile("(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global computed_pattern
computed_pattern = re.compile("computed", flags=re.IGNORECASE)

global comment_pattern
comment_pattern = re.compile('comment.*', flags=re.IGNORECASE)

def load_spectrum_list(msp_file_path):
    """
    Load a spectrum list from a given MSP (Mass Spectral Peak) file.

    :param msp_file_path: The path to the MSP file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.

    Example usage:
    ```
    msp_file_path = "path/to/spectrum.msp"
    spectrum_list = load_spectrum_list(msp_file_path)
    print(spectrum_list)
    ```
    """
    spectrum_list = []
    buffer = []

    total_lines = sum(1 for line in open(msp_file_path, 'r', encoding="UTF-8")) # count the total number of lines in the file

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_lines, unit=" rows", colour="green", desc="\t     reading"): # wrap this with tqdm
            if line.strip() == '':
                if buffer:
                    spectrum_list.append('\n'.join(buffer))
                    buffer = []
            else:
                if not buffer:
                    buffer.append(f"FILENAME: {os.path.basename(msp_file_path)}") # adding filename to spectrum
                buffer.append(line.strip())

    # Add the last spectrum to the list
    if buffer:
        spectrum_list.append('\n'.join(buffer))

    return spectrum_list

def extract_metadata_and_peak_list(spectrum):
    """
    :param spectrum: The spectrum string containing both metadata and peak list
    :return: A tuple containing the extracted metadata and peak list

    Extracts the metadata and peak list from the given spectrum string. The spectrum string should be in a specific format where the metadata appears before the peak list. The method uses
    * regular expressions to find and extract the metadata and peak list.

    The regular expression pattern used to search for the metadata and peak list is:
    "([\s\S]*:.[0-9]*\n)(((-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)"

    The method searches for this pattern in the spectrum string and if a match is found, it extracts the metadata and peak list.

    Example usage:
    spectrum = "Metadata: 123\n1 2\n3 4\n5 6"
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    print(metadata)  # Output: "Metadata: 123"
    print(peak_list)  # Output: "1 2\n3 4\n5 6"
    """
    metadata,peak_list = None,None

    match = re.search(metadata_peak_list_split_pattern, spectrum)

    if match:
        metadata, peak_list = match.group(1), match.group(2)

    return metadata,peak_list

def check_for_metadata_in_comments(metadata_matches):
    """
    :param metadata_matches: A list containing matches of metadata found in comments.
    :return: A new list of metadata matches. Returns False if no metadata matches found.

    This method checks if comment fields exist in each metadata match. If a comment field exists, it further checks if there are sub-fields within the comment field. The sub-fields are extracted
    * using regular expressions and added to the new_metadata_matches list. If no sub-fields are found, the entire match is added to the new_metadata_matches list.

    If no metadata matches are found, the method returns False.

    Example usage:
    metadata_matches = [('comment', 'field1=value1; field2="value2"'), ('metadata', 'field1=value1; field2=value2')]

    result = check_for_metadata_in_comments(metadata_matches)
    print(result)
    # Output: [('comment', 'field1=value1'), ('comment', 'field2="value2"')]

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

def convert_keys(metadata_dict):
    """
    Convert keys in metadata_dict based on the provided keys_dict and keys_list.

    :param metadata_dict: A dictionary containing metadata information.
    :return: A dictionary with converted keys based on the provided keys_dict and keys_list.
    """

    converted = {keys_dict[key]: val for key, val in metadata_dict.items() if
                 key in keys_dict and keys_dict[key] in keys_list}

    converted.update({key: "" for key in keys_list if key not in converted})

    return converted

def metadata_to_dict(metadata):
    """
    Convert metadata string to DataFrame.

    :param metadata: The metadata string to be converted.
    :type metadata: str
    :return: DataFrame containing metadata keys and values.
    :rtype: pandas.DataFrame
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

        metadata_dict = convert_keys(metadata_dict)
        metadata_dict = normalize_values(metadata_dict)

        return metadata_dict

    return metadata_dict

def peak_list_to_df(peak_list, precursormz):
    """
    Converts a peak list string into a numpy array.

    :param peak_list: A string representing a peak list. Each peak is represented by a pair of values, separated by a space or colon. The first value represents the m/z (mass-to-charge ratio
    *) of the peak, and the second value represents the intensity of the peak.
    :return: A numpy array containing the peak data, with two columns for "mz" and "intensity". The "mz" column contains the m/z values, and the "intensity" column contains the corresponding
    * peak intensities.
    """
    peaks_match = re.findall(peak_list_split_pattern, peak_list)

    if peaks_match:
        peaks_match = [(float(i), float(j)) for i, j in peaks_match]

        # Convert list of tuples to numpy array
        peak_array = np.array(peaks_match, dtype=float)

        # Sort the array based on the mz values
        peak_array = peak_array[peak_array[:, 0].argsort()]

        peak_array = apply_filters(peak_array, precursormz)
        return peak_array
    else:
        return np.array([])  # Return an empty numpy array

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    Structure metadata and peak list into formatted DataFrames.

    :param metadata: The metadata to be structured.
    :type metadata: list or tuple
    :param peak_list: The peak list to be structured.
    :type peak_list: list or tuple
    :return: Tuple containing the structured metadata and peak list as DataFrames.
    :rtype: tuple
    """
    if not metadata or not peak_list:
        return {},np.array([])
    else:
        metadata_dict = metadata_to_dict(metadata)
        if not metadata_dict:
            return {},np.array([])
        if "PRECURSORMZ" in metadata_dict:
            if metadata_dict["PRECURSORMZ"]:
                try:
                    peak_list_DF = peak_list_to_df(peak_list,float(metadata_dict["PRECURSORMZ"].replace(",",".")))
                    return metadata_dict, peak_list_DF
                except:
                    return {},np.array([])
            else:
                return {},np.array([])
        else:
            return {},np.array([])

def parse_metadata_and_peak_list(spectrum):
    """
    Parse metadata and peak list from a given spectrum.

    :param spectrum: A spectrum object or data structure that contains both metadata and peak list.
    :return: A tuple containing dataframes for the parsed metadata and peak list.
    """
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    metadata_DF, peak_list_DF = structure_metadata_and_peak_list(metadata, peak_list)

    return metadata_DF, peak_list_DF

def msp_parser(spectrum):
    """
    Parse the given spectrum to extract metadata and peak list.

    :param spectrum: The spectrum data to be parsed.
    :return: The parsed metadata with peak list.
    """
    time.sleep(0.000000001) # Needed to ensure progress bar display update (1ns)

    metadata,peak_list = parse_metadata_and_peak_list(spectrum)

    if metadata == {} or len(peak_list) == 0:
        return None
    else:
        metadata['NUM PEAKS'] = len(peak_list)
        metadata['PEAKS_LIST'] = peak_list
        return metadata

def msp_cleaning_processing(spectrum_list):
    """
    :param spectrum_list: A list of spectrum data to be processed.
    :return: A list containing the results of processing the spectrum data.

    This method takes in a list of spectrum data and processes it using multiple threads for improved performance. The spectrum data is passed to the `msp_parser` function, which carries
    * out the actual processing. The `ThreadPoolExecutor` from the `concurrent.futures` module is used to concurrently execute the `msp_parser` function on each spectrum in the provided
    * list.

    The function returns a list containing the results of processing the spectrum data. Any `None` results are filtered out before returning the final list.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(msp_parser, spectrum_list), total=len(spectrum_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.
