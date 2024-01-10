from msp_to_json import *
from xml_to_json import *
from csv_to_json import *
from tqdm import tqdm
import pandas as pd
import time
import os
import re

def load_spectrum_list_from_msp(msp_file_path):
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
    filename = os.path.basename(msp_file_path)
    spectrum_list = []
    buffer = []

    total_lines = sum(1 for line in open(msp_file_path, 'r', encoding="UTF-8")) # count the total number of lines in the file

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_lines, unit=" rows", colour="green", desc="{:>80}".format(f"loading [{filename}]")): # wrap this with tqdm
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

def concatenate_MSP(msp_list):
    """
    Concatenates a list of MSP files into a single spectrum list.

    :param msp_list: A list of file paths to MSP files.
    :return: A concatenated spectrum list.
    """
    spectrum_list = []

    for files in msp_list:
        spectrum_list.extend(load_spectrum_list_from_msp(files))

    return spectrum_list

def concatenate_xml(xml_list):
    """
    Concatenates the contents of XML files into a single XML string.

    :param xml_list: List of XML file paths to be concatenated.
    :return: List containing the concatenated XML contents of all files.

    Example usage:
        xml_list = ["file1.xml", "file2.xml"]
        result = concatenate_xml(xml_list)
    """
    FINAL_XML = []
    for files in tqdm(xml_list, total=len(xml_list), unit=" spectrums", colour="green", desc="{:>80}".format("concatenate")):
        if files.endswith(".xml"):
            file_name = os.path.basename(files.replace(".xml", ""))
            with open(files, "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            # Add filename to xml
            xml_content = re.sub("</sample-mass>\n",f"</sample-mass>\n  <filename>{file_name}</filename>\n",xml_content)

            FINAL_XML.extend([xml_content])

    return FINAL_XML

def concatenate_csv(csv_list):
    """
    Concatenates a list of CSV files into a single DataFrame.

    :param csv_list: A list of file paths to CSV files.
    :return: A pandas DataFrame containing the concatenated data.

    Example usage:
    ```python
    csv_list = ["file1.csv", "file2.csv", "file3.csv"]
    df = concatenate_csv(csv_list)
    ```
    """
    df_list = [pd.read_csv(file, sep=";", encoding="UTF-8") for file in csv_list]
    df = pd.concat(df_list, ignore_index=True)

    return df





def convert_to_json(input_path):
    """
    Converts JSON and XML files to MSP format.

    :param input_path: The path to the directory containing the JSON and XML files.
    :return: A tuple containing the converted JSON and XML files in MSP format.
    """
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
        time.sleep(0.02)
        print("{:>80}".format("-- CONVERTING MSP TO JSON --"))
        # Concatenate all MSP to a list
        FINAL_MSP = concatenate_MSP(msp_list)
        # Convert all MSP spectrum to JSON spectrum (Multithreaded)
        FINAL_MSP = msp_to_json_processing(FINAL_MSP)

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
        time.sleep(0.02)
        print("{:>80}".format("-- CONVERTING XML TO JSON --"))
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
        time.sleep(0.02)
        print("{:>80}".format("-- CONVERTING CSV TO JSON --"))
        # Concatenate all CSV to a list
        FINAL_CSV = concatenate_csv(csv_list)
        # Convert all CSV spectrum to JSON spectrum (Multithreaded)
        FINAL_CSV = csv_to_json_processing(FINAL_CSV)

    return FINAL_MSP, FINAL_XML, FINAL_CSV
