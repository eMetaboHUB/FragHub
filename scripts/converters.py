import concurrent.futures
import time

from tqdm import tqdm
import json
import os
import re

def concatenate_json(json_path):
    """
    Concatenates multiple JSON files into a single list of dictionaries.

    :param json_path: The directory path containing the JSON files.
    :return: A list of dictionaries containing the concatenated JSON data.
    """
    JSON_LIST = []
    for files in tqdm(os.listdir(json_path), total=len(os.listdir(json_path)), unit=" spectrums", colour="green", desc="\t concatenate"):
        if files.endswith(".json"):
            file_name = os.path.basename(os.path.join(json_path, files)).replace(".json", "")
            with  open(os.path.join(json_path, files), "r", encoding="UTF-8") as f:
                lines = f.readlines()

            data = [json.loads(line) for line in lines]  # returns JSON object as a list of dictionary
            # Add filename to json
            for dicts in data:
                dicts["filename"] = file_name

            JSON_LIST.extend(data)

    return JSON_LIST

def json_to_msp(json_spectrum):
    """
    Convert a JSON spectrum to an MSP (Mass Spectral Peak) string.

    :param json_spectrum: A dictionary representing the JSON spectrum.
                          It should contain the following keys:
                            - "filename" (str): The name of the file.
                            - Any other key-value pairs representing metadata.
                            - "peaks" (list): A list of tuples representing peaks.
                              Each tuple should contain two values:
                                - The peak's m/z value (float)
                                - The peak's intensity (float)

    :return: An MSP string representation of the spectrum.
    """
    SPECTRUM = ""  # Creating empty spectrum string
    SPECTRUM = SPECTRUM + "FILENAME: " + json_spectrum["filename"] + "\n"
    for key, value in json_spectrum.items():
        if key != "peaks":
            SPECTRUM = SPECTRUM + key + ": " + str(value) + "\n"
        else:
            SPECTRUM = SPECTRUM + "num peaks: " + str(len(json_spectrum["peaks"])) + "\n"
            for fragments in json_spectrum["peaks"]:
                SPECTRUM = SPECTRUM + str(fragments[0]) + " " + str(fragments[1]) + "\n"

    return SPECTRUM

def JSON_convert_processing(FINAL_JSON):
    """
    Convert JSON objects to MSP format using multiple threads.

    :param FINAL_JSON: A list of JSON objects to be converted.
    :return: A list containing the results of the conversion.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(json_to_msp, FINAL_JSON), total=len(FINAL_JSON), unit=" spectrums", colour="green", desc="\t  converting"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

def concatenate_xml(xml_path):
    """
    Concatenates the content of XML files in a given directory.

    :param xml_path: The path to the directory containing XML files.
    :return: A list containing the concatenated XML content from all files.
    """
    FINAL_XML = []
    for files in tqdm(os.listdir(xml_path), total=len(os.listdir(xml_path)), unit=" spectrums", colour="green", desc="\t concatenate"):
        if files.endswith(".xml"):
            file_name = os.path.basename(os.path.join(xml_path, files).replace(".xml", ""))
            with open(os.path.join(xml_path, files), "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            # Add filename to xml
            xml_content = re.sub("</sample-mass>\n",f"</sample-mass>\n  <filename>{file_name}</filename>\n",xml_content)

            FINAL_XML.extend([xml_content])

    return FINAL_XML

def xml_to_msp(xml_content):
    """
    Convert XML content to MSP format.

    :param xml_content: The XML content to be converted.
    :return: The converted MSP format.

    """
    specrta_dict = {}

    if re.search("<filename>(.*)</filename>", xml_content):
        specrta_dict["filename"] = re.search("<filename>(.*)</filename>", xml_content).group(1)

    if re.search("<inchikey>(.*)</inchikey>", xml_content):
        specrta_dict["inchikey"] = re.search("<inchikey>(.*)</inchikey>", xml_content).group(1)

    if re.search("<SMILES>(.*)</SMILES>", xml_content):
        specrta_dict["SMILES"] = re.search("<SMILES>(.*)</SMILES>", xml_content).group(1)

    if re.search("<INCHI>(.*)</INCHI>", xml_content):
        specrta_dict["INCHI"] = re.search("<INCHI>(.*)</INCHI>", xml_content).group(1)

    if re.search("<MOLECULAR_FORMULA>(.*)</MOLECULAR_FORMULA>", xml_content):
        specrta_dict["MOLECULAR_FORMULA"] = re.search("<MOLECULAR_FORMULA>(.*)</MOLECULAR_FORMULA>", xml_content).group(1)

    if re.search("<ACC_MASS>(.*)</ACC_MASS>", xml_content):
        specrta_dict["PRECURSORMZ"] = re.search("<ACC_MASS>(.*)</ACC_MASS>", xml_content).group(1)

    if re.search("<instrument-type>(.*)</instrument-type>", xml_content):
        specrta_dict["instrument-type"] = re.search("<instrument-type>(.*)</instrument-type>", xml_content).group(1)

    if re.search("<collision-energy-voltage>(.*)</collision-energy-voltage>", xml_content):
        specrta_dict["collision-energy-voltage"] = re.search("<collision-energy-voltage>(.*)</collision-energy-voltage>", xml_content).group(1)

    if re.search("<ionization-mode>(.*)</ionization-mode>", xml_content):
        specrta_dict["ionization-mode"] = re.search("<ionization-mode>(.*)</ionization-mode>", xml_content).group(1)

    if re.search("<predicted>(.*)</predicted>", xml_content):
        specrta_dict["predicted"] = re.search("<predicted>(.*)</predicted>", xml_content).group(1)

    if re.search("<database-id>(.*)</database-id>", xml_content):
        specrta_dict["database-id"] = re.search("<database-id>(.*)</database-id>", xml_content).group(1)

    # Correcting 0 charge
    if re.search("<ionization-mode>(.*)</ionization-mode>", xml_content):
        if specrta_dict["ionization-mode"][0] == "Positive":
            specrta_dict["charge"] = "1"
        elif specrta_dict["ionization-mode"][0] == "Negative":
            specrta_dict["charge"] = "-1"

    peak_list = [re.findall("<mass-charge>(.*)</mass-charge>",xml_content),re.findall("<intensity>(.*)</intensity>",xml_content)]

    specrta_dict_final = {}

    for key, value in specrta_dict.items():
        specrta_dict_final[key] = None
        if specrta_dict[key] != [] and specrta_dict[key] != None:
            specrta_dict_final[key] = specrta_dict[key]

    # Starting to write information from xml to msp format
    specrta_dict_final["peak_list"] = ""

    for mass_charge, intensity in zip(peak_list[0], peak_list[1]):
        specrta_dict_final["peak_list"] = specrta_dict_final["peak_list"] + mass_charge + " " + intensity + "\n"

    SPECTRUM = ""

    for key, value in specrta_dict_final.items():
        if key == "peak_list":
            SPECTRUM = SPECTRUM + "NUM PEAKS: " + str(len(peak_list[0])) + "\n"
            SPECTRUM = SPECTRUM + str(specrta_dict_final[key])
        else:
            SPECTRUM = SPECTRUM + key + ": " + str(specrta_dict_final[key]) + "\n"

    return SPECTRUM

def XML_convert_processing(FINAL_XML):
    """
    :param FINAL_XML: List of XML elements to be converted to MSP format.
    :return: List of converted XML elements in MSP format.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(xml_to_msp, FINAL_XML), total=len(FINAL_XML), unit=" spectrums", colour="green", desc="\t  converting"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


def convert_to_msp(input_path):
    """
    Convert JSON and XML files to MSP format.

    :param input_path: The path to the directory containing JSON and XML files.
    :return: A tuple containing the final JSON and XML converted to MSP format.
    """
    # JSON
    FINAL_JSON = []
    json_to_do = False
    json_path = os.path.join(input_path,"JSON")
    # check if there is a json file into the directory
    for files in os.listdir(json_path):
        if files.endswith(".json"):
            json_to_do = True
    if json_to_do == True:
        time.sleep(0.02)
        print("-- CONVERTING JSON TO MSP --")
        # Concatenate all JSON to a list
        FINAL_JSON = concatenate_json(json_path)
        # Convert all JSON spectrum to MSP spectrum (Multithreaded)
        FINAL_JSON = JSON_convert_processing(FINAL_JSON)


    # XML
    FINAL_XML = []
    xml_to_do = False
    xml_path = os.path.join(input_path,"XML")
    # check if there is a xml file into the directory
    for files in os.listdir(xml_path):
        if files.endswith(".xml"):
            xml_to_do = True
    if xml_to_do == True:
        time.sleep(0.02)
        print("-- CONVERTING XML TO MSP --")
        # Concatenate all XML to a list
        FINAL_XML = concatenate_xml(xml_path)
        # Convert all XML spectrum to MSP spectrum (Multithreaded)
        FINAL_XML = XML_convert_processing(FINAL_XML)

    return FINAL_JSON,FINAL_XML