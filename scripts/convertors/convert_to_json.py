from .json_to_json import *
from .msp_to_json import *
from .xml_to_json import *
from .csv_to_json import *
from .mgf_to_json import *
from .loaders import *
from tqdm import tqdm
import pandas as pd
import time
import json
import os
import re

def concatenate_MSP(msp_list):
    """
    Concatenates multiple MSP files into a single spectrum list.

    :param msp_list: A list of paths to MSP files.
    :return: A list of spectra loaded from the MSP files.
    """
    spectrum_list = []

    for files in msp_list:
        spectrum_list.extend(load_spectrum_list_from_msp(files))

    return spectrum_list

def concatenate_xml(xml_list):
    """
    Concatenates XML files and adds filename tags to each file's XML content.

    :param xml_list: List of XML file paths.
    :return: List of concatenated XML contents with added filename tags.
    """
    FINAL_XML = []
    for files in tqdm(xml_list, total=len(xml_list), unit=" spectrums", colour="green", desc="{:>70}".format("concatenate")):
        if files.endswith(".xml"):
            file_name = os.path.basename(files.replace(".xml", ""))
            with open(files, "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            # Add filename to xml
            xml_content = re.sub("</id>\n",f"</id>\n  <filename>{file_name}</filename>\n",xml_content)

            FINAL_XML.extend([xml_content])

    return FINAL_XML

def concatenate_csv(csv_list):
    """
    Concatenates multiple CSV files into a single DataFrame.

    :param csv_list: List of paths to CSV files to be concatenated.
    :type csv_list: list[str]
    :return: Concatenated DataFrame.
    :rtype: pandas.DataFrame
    """
    df_list = []

    for file in csv_list:
        df = pd.read_csv(file, sep=";", quotechar='"', encoding="UTF-8", dtype=str)
        df['filename'] = os.path.basename(file)

        df.columns = df.columns.str.lower()

        df = df.astype(str)

        df_list.append(df)

    df = pd.concat(df_list, ignore_index=True)

    return df

def concatenate_JSON(json_list):
    """
    Concatenates a list of JSON files into a single spectrum list.

    :param json_list: A list of JSON file paths.
    :type json_list: list[str]
    :return: The concatenated spectrum list.
    :rtype: list[dict]
    """
    spectrum_list = []

    for files in json_list:
        spectrum_list.extend(load_spectrum_list_json(files))

    return spectrum_list

def concatenate_MGF(mgf_list):
    """
    Concatenates multiple MGF files into a single spectrum list.

    :param mgf_list: A list of paths to MGF files.
    :type mgf_list: list[str]
    :return: The concatenated spectrum list.
    :rtype: list[str]
    """
    spectrum_list = []

    for files in mgf_list:
        spectrum_list.extend(load_spectrum_list_from_mgf(files))

    return spectrum_list

def convert_to_json(input_path):
    """
    :param input_path: The path of the input directory where the files are located.
    :return: A tuple containing the converted data from each file type to JSON format. The order of the elements in the tuple is as follows:
        - FINAL_MSP: The converted data from MSP files to JSON format.
        - FINAL_XML: The converted data from XML files to JSON format.
        - FINAL_CSV: The converted data from CSV files to JSON format.
        - FINAL_JSON: The converted data from JSON files to JSON format (not actually converted, just collected).
        - FINAL_MGF: The converted data from MGF files to JSON format.
    """
    # JSON
    FINAL_JSON = []
    json_list = []
    json_to_do = False
    json_path = os.path.join(input_path, "JSON")
    # check if there is a json file into the directory
    for root, dirs, files in os.walk(json_path):
        for file in files:
            if file.endswith(".json"):
                json_path = os.path.join(root, file)  # Full path to the file
                json_list.append(json_path)
                json_to_do = True
    if json_to_do == True:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING JSON TO JSON --"))
        time.sleep(0.01)
        # Concatenate all JSON to a list
        FINAL_JSON = concatenate_JSON(json_list)
        # Convert all bad structured JSON to pretty structured JSON (Multithreaded)
        FINAL_JSON = json_to_json_processing(FINAL_JSON)

    # MSP
    FINAL_MSP = []
    msp_list = []
    msp_to_do = False
    msp_path = os.path.join(input_path,"MSP")
    # check if there is a json file into the directory
    for root, dirs, files in os.walk(msp_path):
        for file in files:
            if file.endswith(".msp"):
                msp_path = os.path.join(root, file)  # Full path to the file
                msp_list.append(msp_path)
                msp_to_do = True
    if msp_to_do == True:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING MSP TO JSON --"))
        time.sleep(0.01)
        # Concatenate all MSP to a list
        FINAL_MSP = concatenate_MSP(msp_list)
        # Convert all MSP spectrum to JSON spectrum (Multithreaded)
        FINAL_MSP = msp_to_json_processing(FINAL_MSP)

    # MGF
    FINAL_MGF = []
    mgf_list = []
    mgf_to_do = False
    mgf_path = os.path.join(input_path, "MGF")
    # check if there is a json file into the directory
    for root, dirs, files in os.walk(mgf_path):
        for file in files:
            if file.endswith(".mgf"):
                mgf_path = os.path.join(root, file)  # Full path to the file
                mgf_list.append(mgf_path)
                mgf_to_do = True
    if mgf_to_do == True:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING MGF TO JSON --"))
        time.sleep(0.01)
        # Concatenate all MGF to a list
        FINAL_MGF = concatenate_MGF(mgf_list)
        # Convert all MGF spectrum to JSON spectrum (Multithreaded)
        FINAL_MGF = mgf_to_json_processing(FINAL_MGF)

    # XML
    FINAL_XML = []
    xml_list = []
    xml_to_do = False
    xml_path = os.path.join(input_path, "XML")
    # check if there is a xml file into the directory
    for root, dirs, files in os.walk(xml_path):
        for file in files:
            if file.endswith(".xml"):
                xml_path = os.path.join(root, file)  # Full path to the file
                xml_list.append(xml_path)
                xml_to_do = True
    if xml_to_do == True:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING XML TO JSON --"))
        time.sleep(0.01)
        # Concatenate all XML to a list
        FINAL_XML = concatenate_xml(xml_list)
        # Convert all XML spectrum to XML spectrum (Multithreaded)
        FINAL_XML = xml_to_json_processing(FINAL_XML)

    # CSV
    FINAL_CSV = []
    csv_list= []
    csv_to_do = False
    csv_path = os.path.join(input_path, "CSV")
    # check if there is a csv file into the directory
    for root, dirs, files in os.walk(csv_path):
        for file in files:
            if file.endswith(".csv"):
                csv_path = os.path.join(root, file)  # Full path to the file
                csv_list.append(csv_path)
                csv_to_do = True
    if csv_to_do == True:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING CSV TO JSON --"))
        time.sleep(0.01)
        # Concatenate all CSV to a list
        FINAL_CSV = concatenate_csv(csv_list)
        # Convert all CSV spectrum to JSON spectrum (Multithreaded)
        FINAL_CSV = csv_to_json_processing(FINAL_CSV)

    return FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF
