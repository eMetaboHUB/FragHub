import concurrent.futures
from tqdm import tqdm
from loaders import *
import hashlib
import ijson
import json
import os
import re

global inchikey_update_pattern
inchikey_update_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])", flags=re.IGNORECASE)

global peak_list_update_pattern
peak_list_update_pattern = re.compile(r"\"peaks\":([\S\s]*?)\}", flags=re.IGNORECASE)

global peak_list_split_update_pattern
peak_list_split_update_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def hash_spectrum_data(spectrum_data):
    """
    Calculate the SHA256 hash of spectrum data.

    :param spectrum_data: The spectrum data to hash.
    :type spectrum_data: str
    :return: The SHA256 hash of spectrum data.
    :rtype: str
    """
    spectrum_string = str(spectrum_data)

    inchikey = re.search(inchikey_update_pattern, spectrum_string)
    peak_list = re.search(peak_list_update_pattern, spectrum_string)

    if inchikey:
        inchikey = inchikey.group(1)

    if inchikey and peak_list:
        peak_list = peak_list.group(1)
        spectrum_string = inchikey+"\n"+peak_list

    else:
        spectrum_string = str(spectrum_data)

    # Créer un objet sha256
    sha256 = hashlib.sha256()

    # Fournir les données de spectre à sha256
    sha256.update(spectrum_string.encode('utf-8'))

    # Retourner le hash sha256 en hex
    return sha256.hexdigest()

def genrate_fraghubid(spectrum):
    """
    Generate FragHubID for a given spectrum.

    :param spectrum: The spectrum data.
    :return: The spectrum data with Fragment Hub ID.
    """
    fraghubid = str(hash_spectrum_data(spectrum))
    new_spectrum = {"FRAGHUBID": fraghubid}
    new_spectrum.update(spectrum)

    return new_spectrum

def genrate_fraghubid_processing(spectrum_list, files):
    """
    Perform parallel processing of the given spectrum list and generate fraghubid for each spectrum.

    :param spectrum_list: A list of spectra.
    :return: A list of fraghubids generated for each spectrum.
    """
    filename = os.path.basename(files)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(genrate_fraghubid, spectrum_list), total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>80}".format(f"generating FragHubID on [{filename}]")))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


def check_fraghubid_already_done(json_file_path):
    """
    Check if a FragHubID is already done in a given file.

    :param json_file_path: The file path of the JSON file to check.
    :type json_file_path: str
    :return: True if a FRAGHUBID is found, False otherwise.
    :rtype: bool
    """
    with open(json_file_path, 'r') as file:
        # ijson items returns a generator yielding items in a json file
        objects = ijson.items(file, '')

        first_dict = next(objects)

        # check if 'FRAGHUBID' is in the first dictionary
        if 'FRAGHUBID' in first_dict:
            return True

    # 'FRAGHUBID' was not found or file was empty
    return False

def process_converted_after(spectrum_list, mode):
    """
    :param spectrum_list: List of spectra to be processed
    :param mode: Mode in which the spectra are converted (MSP, XML, CSV)
    :return: None

    Processes the converted spectra and writes them to a JSON file based on the specified mode.

    The `spectrum_list` parameter should contain a list of spectra to be processed. The `mode` parameter specifies the format in which the spectra are converted, and can be one of the following
    *: "MSP", "XML", or "CSV".

    The method first determines the file path and filename based on the mode. The file path is set to the absolute path of the corresponding converted JSON file, and the filename is extracted
    * from the file path.

    Next, the method calls the `genrate_fraghubid_processing` function to generate additional processing on the spectrum list. The `spectrum_list` and `filename` parameters are passed to
    * this function, which modifies and returns the spectrum list.

    If the spectrum list is not empty, the method opens the file specified by the file path in write mode and writes the spectra to the file in JSON format. It iterates over the spectrum
    * list and writes each spectrum object as a JSON object to the file. The JSON objects are separated by commas. Finally, the method writes closing brackets to end the JSON list.

    Note: The method uses the `tqdm` package to display a progress bar while writing the spectra to the file.
    """
    if mode == "MSP":
        file_path = os.path.abspath("../INPUT/JSON/MSP_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "XML":
        file_path = os.path.abspath("../INPUT/JSON/XML_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "CSV":
        file_path = os.path.abspath("../INPUT/JSON/CSV_converted.json")
        filename = os.path.basename(file_path)

    spectrum_list = genrate_fraghubid_processing(spectrum_list, filename)

    if spectrum_list:
        with open(file_path, "w", encoding="UTF-8") as buffer:
            buffer.write("[")  # begin the JSON list
            for i in tqdm(range(len(spectrum_list)), total=len(spectrum_list), unit=" row", colour="green", desc="{:>80}".format(f"writting {filename}")):
                # write a comma before every object except the first one
                if i != 0:
                    buffer.write(",")
                json.dump(spectrum_list[i], buffer, ensure_ascii=False)
            buffer.write("]")  # end the JSON lists

def generate_fraghub_id(json_directory_path, FINAL_MSP, FINAL_XML, FINAL_CSV):
    """
    :param msp_directory_path: The path to the MSP directory.
    :return: None

    This method generates a FragHub ID for each MSP file in the given directory. It loops through all the files in the directory and checks if they have a .msp extension. For each .msp file
    *, it checks if a FragHub ID has already been generated by calling the `check_fraghubid_already_done` method. If a FragHub ID hasn't been generated, it loads the spectrum list from the
    * MSP file using the `load_spectrum_list` method. Then, it processes the spectrum list to generate the FragHub ID using the `generate_fraghubid_processing` method.

    Finally, it opens the MSP file in write mode and writes the modified spectrum list with the generated FragHub ID into the file.
    """
    json_path_list = []

    for root, dirs, files in os.walk(json_directory_path):
        for file in files:
            if file.endswith(".json"):
                json_path_list.append(os.path.join(root, file))

    for files in json_path_list:
        if files.endswith(".json"):
            filename = os.path.basename(files)
            if not check_fraghubid_already_done(str(files)):
                spectrum_list = load_spectrum_list_json(files)
                spectrum_list = genrate_fraghubid_processing(spectrum_list, files)

                if spectrum_list:
                    with open(files, "w", encoding="UTF-8") as buffer:
                        buffer.write("[")  # begin the JSON list
                        for i in tqdm(range(len(spectrum_list)), total=len(spectrum_list), unit=" row", colour="green", desc="{:>80}".format(f"writting {filename}")):
                            # write a comma before every object except the first one
                            if i != 0:
                                buffer.write(",")
                            json.dump(spectrum_list[i], buffer, ensure_ascii=False)
                        buffer.write("]")  # end the JSON list

    if FINAL_MSP:
        process_converted_after(FINAL_MSP, "MSP")
    if FINAL_XML:
        process_converted_after(FINAL_XML, "XML")
    if FINAL_CSV:
        process_converted_after(FINAL_CSV, "CSV")