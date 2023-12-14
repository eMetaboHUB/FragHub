import concurrent.futures
from tqdm import tqdm
import pandas as pd
import time
import re


def load_spectrum_list(msp_file_path):
    """
    Load spectrum list from a given MSP file.

    :param msp_file_path: The file path of the MSP file to be loaded.
    :return: A list of spectra. Each spectrum is represented as a string.
    """
    spectrum_list = []
    buffer = []

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in file:
            if line.strip() == '':
                if buffer:  # Only add the buffer to the list if it's not empty
                    spectrum_list.append('\n'.join(buffer))
                    buffer = []
            else:
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

    match = re.search("([\s\S]*:.[0-9]*\n)(((-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)", spectrum)

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
        if re.search("comment.*", match[0], flags=re.IGNORECASE):
            if "=" in match[1]:
                sub_fields_matches = re.findall('(\S+)=\"([^\"]*)\"|\"(\w+)=([^\"]*)\"|\"([^\"]*)=([^\"]*)\"|(\S+)=(\d+(?:\.\d*)?)|(\S+?)=(.*?)(;|\n|$)', match[1])
                if sub_fields_matches:
                    for sub_fields_match in sub_fields_matches:
                        non_empty_tuple = tuple(group for group in sub_fields_match if group)
                        if non_empty_tuple:
                            new_metadata_matches.append(non_empty_tuple)
                else:
                    return False
            else:
                new_metadata_matches.append(match)
        else:
            new_metadata_matches.append(match)

    return new_metadata_matches if new_metadata_matches else False

def metadata_to_df(metadata):
    """
    Convert metadata string to DataFrame.

    :param metadata: The metadata string to be converted.
    :type metadata: str
    :return: DataFrame containing metadata keys and values.
    :rtype: pandas.DataFrame
    """
    metadata_dict = {}

    metadata_matches = re.findall("([\s\S]*?):([\s\S]*?)(\n|$)",metadata)

    if metadata_matches:
        temp = check_for_metadata_in_comments(metadata_matches)
        if temp != False:
            metadata_matches = temp

        for match in metadata_matches:
            metadata_dict[re.sub(r'^[\W_]+|[\W_]+$', '', match[0]).lower()] = [match[1]]

        df = pd.DataFrame.from_dict(metadata_dict)

        return df
    else:
        return pd.DataFrame()

def peak_list_to_df(peak_list):
    """
    Converts a peak list string into a pandas DataFrame.

    :param peak_list: A string representing a peak list. Each peak is represented by a pair of values, separated by a space or colon. The first value represents the m/z (mass-to-charge ratio
    *) of the peak, and the second value represents the intensity of the peak.
    :return: A pandas DataFrame containing the peak data, with two columns named "mz" and "intensity". The "mz" column contains the m/z values, and the "intensity" column contains the corresponding
    * peak intensities.
    """
    peaks_match = re.findall("(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)", peak_list)
    if peaks_match:
        peaks_match = [(float(i), float(j)) for i, j in peaks_match]
        peak_list_DF = pd.DataFrame(peaks_match, columns=["mz", "intensity"])
        return peak_list_DF
    else:
        return pd.DataFrame()

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
    if metadata == None or peak_list == None:
        return pd.DataFrame(),pd.DataFrame()
    else:
        metadata_DF = metadata_to_df(metadata)
        peak_list_DF = peak_list_to_df(peak_list)
        return metadata_DF, peak_list_DF

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

    if metadata.empty or peak_list.empty:
        return None
    else:
        metadata['peak_list'] = [peak_list.copy()]
        return metadata

def msp_parsing_processing(spectrum_list):
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
