from .HMDB_completion import *
import concurrent.futures
from tqdm import tqdm
import re

global metadata_pattern
metadata_pattern = re.compile(r"^  <(.*?)>(.*?)(?:</.*?>)\n")

global peak_list_pattern
peak_list_pattern = re.compile(r"^      <mass-charge>(.*?)</mass-charge>\n      <intensity>(.*?)</intensity>")

global hmdb_pattern
hmdb_pattern = re.compile(r"HMDB.*")

def structure_metadata_and_peak_list(metadata, peak_list):
    """
    :param metadata: the metadata of the spectrum
    :param peak_list: the list of peaks in the spectrum
    :return: the structured spectrum JSON object
    """
    spectrum_json = {key.lower(): value for (key, value) in metadata}
    spectrum_json["peaks"] = [[i, j] for i, j in peak_list]

    return spectrum_json

def extract_metadata_and_peak_list(spectrum):
    """
    Extracts metadata and peak list from a given spectrum.

    :param spectrum: The spectrum to extract metadata and peak list from.
    :type spectrum: str
    :return: A tuple containing the extracted metadata and peak list.
    :rtype: tuple
    """
    metadata = re.findall(metadata_pattern, spectrum)
    peak_list = re.findall(peak_list_pattern, spectrum)

    return metadata, peak_list

def complete_HMDB(spectrum_json):
    """
    Completes the HMDB information in the spectrum JSON.

    :param spectrum_json: The spectrum JSON object to complete.
    :type spectrum_json: dict
    :return: The completed spectrum JSON object.
    :rtype: dict
    """
    if re.search(hmdb_pattern, spectrum_json["filename"]):
        spectrum_json = complete_HMDB_spectrum(spectrum_json)

    return spectrum_json

def xml_to_json(spectrum):
    """
    Converts XML formatted spectrum data to JSON format.

    :param spectrum: XML formatted spectrum data.
    :return: JSON formatted spectrum data.
    """
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)
    spectrum_json = structure_metadata_and_peak_list(metadata, peak_list)
    spectrum_json = complete_HMDB(spectrum_json)

    return spectrum_json

def xml_to_json_processing(FINAL_XML):
    """
    Converts XML spectrums to JSON format.

    :param FINAL_XML: The list of XML spectrums to be converted.
    :return: The list of JSON spectrums after conversion.
    """
    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(FINAL_XML), unit=" spectrums", colour="green", desc="{:>80}".format("converting XML spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(FINAL_XML), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_XML[i:i + chunk_size]
            results = list(executor.map(xml_to_json, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final