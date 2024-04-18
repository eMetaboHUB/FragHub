from .keys_convertor import *
from .HMDB_completion import *
import concurrent.futures
from tqdm import tqdm
import re

global metadata_pattern
metadata_pattern = re.compile(r"^  <(.*?)>(.*?)(?:</.*?>)\n", flags=re.MULTILINE)

global peak_list_pattern
peak_list_pattern = re.compile(r"^      <mass-charge>(.*?)</mass-charge>\n      <intensity>(.*?)</intensity>", flags=re.MULTILINE)

global hmdb_pattern
hmdb_pattern = re.compile(r"HMDB.*")

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    This function is responsible for structuring metadata and peak list into a spectrum JSON object.

    :param metadata: the metadata of the spectrum
    This metadata is a dictionary where keys are the metadata attributes and values are the corresponding metadata values.

    :param peak_list: the list of peaks in the spectrum
    This peak list is a list of tuples/lists where each tuple/list contains the m/z ratio and intensity of the peak.

    :return: the structured spectrum JSON object
    This is a JSON object derived from the metadata and peak list input which can be used for downstream tasks such as storage, transmission and computation.
    """
    # Convert all keys in the metadata dictionary to lowercase.
    spectrum_json = {key.lower(): value for (key, value) in metadata.items()}

    # Convert the peak list into a list of lists where each sub-list contains the m/z ratio and intensity of the peak.
    spectrum_json["peaks"] = [[float(i), float(j)] for i, j in peak_list]

    # Return the structured spectrum JSON object.
    return spectrum_json

def extract_metadata_and_peak_list(spectrum):
    """
    This function takes a spectrum as input and extracts metadata and peak list from it.
    Here, spectrum should be a string. The function uses regular expression patterns (metadata_pattern and
    peak_list_pattern) to find all the matches in the given spectrum and returns a tuple that contains
    the extracted metadata and peak list.

    :param spectrum: The spectrum to extract metadata and peak list from.
    :type spectrum: str
    :return: A tuple containing the extracted metadata and peak list.
    :rtype: tuple
    """
    # Using the pre-defined 'metadata_pattern' regex to extract all metadata from the spectrum
    metadata = re.findall(metadata_pattern, spectrum)

    # Using the pre-defined 'peak_list_pattern' regex to extract all peak list from the spectrum
    peak_list = re.findall(peak_list_pattern, spectrum)

    # Returning a tuple of the extracted metadata and peak list
    return metadata, peak_list

def complete_HMDB(spectrum_json):
    """
    Completes the HMDB information in the spectrum JSON.
    :param spectrum_json: The spectrum JSON object to complete.
    :type spectrum_json: dict
    :return: The completed spectrum JSON object.
    :rtype: dict
    """
    # If the filename of the given spectrum_json object matches the HMDB pattern, call the function 'complete_HMDB_spectrum'
    # to complete the HMDB information. If not, the spectrum_json object will remain unmodified.
    if re.search(hmdb_pattern, spectrum_json["filename"]):
        spectrum_json = complete_HMDB_spectrum(spectrum_json)

    # returning the modified or unmodified spectrum_json object
    return spectrum_json

def xml_to_json(spectrum):
    """
    This function is responsible for converting the spectrum data from XML format to JSON.

    Steps involved:
    1. Separation of metadata and peak list from the spectrum data.
    2. Structuring the metadata and peak list in JSON format
    3. If the file matches the HMDB pattern, then completing the HMDB spectrum with additional information.
    4. Converting the metadata keys to match with the keys present in keys_dict.
       If keys from keys_list are not present after conversion, they are added with an empty value.

    :param spectrum: A piece of information on the spectrum stored in XML format.
    :return: spectrum information in JSON format
    """

    # Splitting the XML spectrum to metadata and peak list
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)

    # Structuring the metadata and peak list to form JSON formatted spectrum data
    spectrum_json = structure_metadata_and_peak_list(metadata, peak_list)

    # If the file is of HMDB pattern, then completing the additional HMDB information in it.
    spectrum_json = complete_HMDB(spectrum_json)

    # Transforming the metadata keys to match with keys dictionary and adding missing keys from keys list.
    spectrum_json = convert_keys(spectrum_json)

    # Returning the JSON formatted spectrum data
    return spectrum_json

def xml_to_json_processing(FINAL_XML):
    """
    Converts XML spectrums to JSON format.
    :param FINAL_XML: The list of XML spectrums to be converted.
    :return: The list of JSON spectrums after conversion.
    """
    chunk_size = 5000  # Define chunk size for concurrent execution
    final = []  # Initiate an empty list to collect JSON results

    # Initiate a progress bar to observe the process of conversion
    progress_bar = tqdm(total=len(FINAL_XML), unit=" spectrums", colour="green", desc="{:>70}".format("converting XML spectrums"))

    # the main loop that divides the XML list into chunks
    for i in range(0, len(FINAL_XML), chunk_size):
        # Use concurrent futures to execute processing in parallel
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_XML[i:i + chunk_size]  # Form a chunk from the list of XMLs
            results = list(executor.map(xml_to_json, chunk))  # Apply xml_to_json function to chunks in parallel
            progress_bar.update(len(chunk))  # Update the progress bar with the chunk size

        # Extend the final list with non-None results from the processed chunk
        final.extend([res for res in results if res is not None])

    progress_bar.close()  # Close the progress bar after the processing is complete
    return final  # Return the list of converted JSONs
