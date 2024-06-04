from splash import Spectrum, SpectrumType, Splash
from convertors.loaders import *
import concurrent.futures
from tqdm import tqdm
import hashlib
import ijson
import json
import sys
import os

def hash_spectrum_data(spectrum_data):
    """
    Calculate the SHA256 hash of spectrum data.
    :param spectrum_data: The spectrum data to hash.
    :type spectrum_data: str
    :return: The SHA256 hash of spectrum data.
    :rtype: str
    """
    # Convert spectrum data to string

    # Search for inchikey_update_pattern in spectrum data
    print(spectrum_data["PEAKS_LIST"])
    peak_list = str(spectrum_data["PEAKS_LIST"])

    # Check if both inchikey and peak_list exist
    if peak_list:
        # Combine inchikey and peak list into one string with a newline separator
        spectrum_string =  "\n" + peak_list

    # Create a new sha256 hash object
    sha256 = hashlib.sha256()

    # Update the sha256 object with the spectrum data
    sha256.update(spectrum_string.encode('utf-8'))

    # Return the hexadecimal representation of the sha256 hash
    return sha256.hexdigest()

def genrate_fraghubid(spectrum):
    """
    This function generates FragHubID for a given spectrum by hashing the spectrum data.

    :param spectrum: The spectrum data.
    :return: The spectrum data with Fragment Hub ID.
    """
    # Hash the spectrum data and convert the resultant hash into a string.
    # fraghubid holds the hashed id of the spectrum data
    fraghubid = str(hash_spectrum_data(spectrum))

    # Return None if fraghubid is empty i.e., no FragHubID could be generated.
    if not fraghubid:
        return None

    # Add the generated fraghubid to the spectrum data dictionary.
    # This line is updating the spectrum dictionary with the newly generated fraghubid
    spectrum["FRAGHUBID"] = fraghubid

    # Return the spectrum data with the newly added FragHubID
    # This will be a dictionary containing the spectrum data, updated with a new FragHubID.
    return spectrum

def generate_fraghubid_processing(spectrum_list, files):
    """
    Perform parallel processing of the given spectrum list and generate fraghubid for each spectrum.
    :param spectrum_list: A list of spectra.
    :return: A list of fraghubids generated for each spectrum.
    """
    filename = os.path.basename(files)  # Extract the base name of the file path
    chunk_size = 5000  # Set the size of chunks
    final = []  # Create an empty list to store final results
    # Create a progress bar with the total length of spectrum_list and a relevant description
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>70}".format(f"generating FragHubID on [{filename}]"))
    # Dividing the spectrum list into chunks
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]  # create a chunk
            # Use a executor's map function to apply 'genrate_fraghubid' function to each spectrum in the chunk, and convert it to a list
            results = list(executor.map(genrate_fraghubid, chunk))
            progress_bar.update(len(chunk))  # update the progress bar by the size of the chunk processed
        final.extend([res for res in results if res is not None])  # Extend the final results list with results that are not None
    progress_bar.close()  # close the progress bar
    return final  # return the final results

# process_converted_after function processes the converted spectrum list.
# Based on the mode, it identifies the appropriate file path for the converted file.
def process_converted_after(spectrum_list, mode):
    """
    Process converted spectrum list after conversion.
    :param spectrum_list: List of converted spectra.
    :param mode: Conversion mode. Can be "MSP", "XML", "CSV", or "JSON".
    :return: Processed spectrum list.
    """
    # Assign relevant file path based on conversion mode
    if mode == "MSP":
        file_path = os.path.abspath("../INPUT/CONVERTED/MSP_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "XML":
        file_path = os.path.abspath("../INPUT/CONVERTED/XML_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "CSV":
        file_path = os.path.abspath("../INPUT/CONVERTED/CSV_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "JSON":
        file_path = os.path.abspath("../INPUT/CONVERTED/JSON_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "MGF":
        file_path = os.path.abspath("../INPUT/CONVERTED/MGF_converted.json")
        filename = os.path.basename(file_path)

    # Generate FragHubID for the spectrum list
    spectrum_list = generate_fraghubid_processing(spectrum_list, filename)

    # Return the processed spectrum list
    return spectrum_list

def generate_fraghub_id(FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF):
    """
    Process the converted files and generate a FragHub ID for each file type.
    - FINAL_MSP: The path to the converted MSP file.
    - FINAL_XML: The path to the converted XML file.
    - FINAL_CSV: The path to the converted CSV file.
    - FINAL_JSON: The path to the converted JSON file.
    - FINAL_MGF: The path to the converted MGF file.

    - Return: A tuple containing the FragHub IDs for each file type (MSP, XML, CSV, JSON, MGF).
    """

    # Check if the MSP file path is valid, if it is, process the file to generate fraghub id
    if FINAL_MSP:
        FINAL_MSP = process_converted_after(FINAL_MSP, "MSP")

    # Similar check for XML file path
    if FINAL_XML:
        FINAL_XML = process_converted_after(FINAL_XML, "XML")

    # Similar check for CSV file path
    if FINAL_CSV:
        FINAL_CSV = process_converted_after(FINAL_CSV, "CSV")

    # Similar check for JSON file path
    if FINAL_JSON:
        FINAL_JSON = process_converted_after(FINAL_JSON, "JSON")

    # Similar check for MGF file path
    if FINAL_MGF:
        FINAL_MGF = process_converted_after(FINAL_MGF, "MGF")

    # Finally return a tuple containing fraghub ids for each file type
    return FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF
