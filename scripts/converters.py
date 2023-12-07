import concurrent.futures
from tqdm import tqdm
import pandas as pd
import numpy as np
import chardet
import time
import json
import os
import re

def concatenate_json(json_path):
    """
    :param json_path: The path to the directory containing the JSON files.
    :return: A list of JSON objects, with each object representing the contents of a JSON file in the given directory.
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
    Converts a JSON spectrum to an MSP format.

    :param json_spectrum: The JSON spectrum to convert.
    :type json_spectrum: dict
    :return: The spectrum in MSP format.
    :rtype: str
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
    Process JSON data and convert it using multiple worker threads.

    :param FINAL_JSON: The JSON data to be processed and converted.
    :type FINAL_JSON: list
    :return: The list of different worker executions.
    :rtype: list
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(json_to_msp, FINAL_JSON), total=len(FINAL_JSON), unit=" spectrums", colour="green", desc="\t  converting"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

def concatenate_xml(xml_path):
    """
    :param xml_path: The directory path where the XML files are located.
    :return: The concatenated XML content in a list.
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

def detect_encoding(file_path):
    """
    Detects the encoding of a file.

    :param file_path: The path to the file.
    :return: The encoding of the file.
    """
    with open(file_path, 'rb') as f:
        result = chardet.detect(f.read())
    return result['encoding']

def try_separators(file_path, encoding):
    """
    Tries different separators to parse a CSV file.

    :param file_path: The path to the input CSV file.
    :param encoding: The encoding of the CSV file.
    :return: A tuple containing the parsed DataFrame and the successful separator.
    :raises ValueError: If the input file could not be parsed with any of the attempted separators.
    """
    separators = [',', ';', '\t', '|']
    for sep in separators:
        try:
            df = pd.read_csv(file_path, sep=sep, encoding=encoding)
            return df
        except pd.errors.ParserError:
            continue
    raise ValueError(f'Input file at {file_path} could not be parsed with the attempted separators.')

def concatenate_csv(csv_path):
    """
    Concatenates all CSV files in a given directory.

    :param csv_path: The path to the directory containing the CSV files.
    :return: A list of dataframes representing the concatenated CSV files.
    """
    FINAL_CSV = []
    for files in tqdm(os.listdir(csv_path), total=len(os.listdir(csv_path)), unit=" spectrums", colour="green", desc="\t concatenate"):
        if files.endswith(".csv"):
            file_path = os.path.join(csv_path, files)
            file_name = os.path.basename(file_path.replace(".csv", ""))
            encoding = detect_encoding(file_path)
            csv_df = try_separators(file_path, encoding)

            # Add filename column to dataframe
            csv_df.insert(0, 'filename', file_name)

            FINAL_CSV.extend([csv_df])

    return FINAL_CSV

def xml_to_msp(xml_content):
    """
    Converts XML content to MSP format.

    :param xml_content: XML content to be converted.
    :return: Converted content in MSP format.
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
    Process XML conversion using multithreading.

    :param FINAL_XML: The list of XML files to be converted.
    :return: The list of successfully converted XML files.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(xml_to_msp, FINAL_XML), total=len(FINAL_XML), unit=" spectrums", colour="green", desc="\t  converting"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

def csv_to_msp(dataframe):
    """
    Converts a DataFrame to a formatted string with the following format for each row:
    {column_name}: {value}
    The rows are concatenated with newline characters.

    :param dataframe: A pandas DataFrame object.
    :return: A string with formatted rows for each entry in the DataFrame.
    """

    def format_row(row):
        """
        :param row: The row of a DataFrame to format.
        :return: The formatted row as a string, with each column and value represented as `<column>: <value>` or just `<value>` if the value is a peak list.

        """
        return '\n'.join([f'{col}: {val}' if col != 'PEAK_LIST' else f'{val}' for col, val in zip(row.index, row.values)])

    def is_peak(val):
        """
        :param val: The value to be checked if it represents a peak. The value can be a string, integer, or float.
        :return: Returns True if the value represents a peak, otherwise False. A value represents a peak if it matches the pattern '(\d+\.?\d*)\s+(\d+\.?\d*)', where the first group represents
        * the peak's x-coordinate and the second group represents the peak's y-coordinate.
        """
        return bool(re.search(r'(\d+\.?\d*)\s+(\d+\.?\d*)', str(val)))

    peak_mask = dataframe.applymap(is_peak)
    # Identify columns where at least one element matches the regex
    peak_columns = peak_mask.all()

    for col in peak_columns.index:
        if peak_columns[col]:
            dataframe.rename(columns={col: 'PEAK_LIST'}, inplace=True)

    # apply format_row to each row of the DataFrame
    formatted_rows = dataframe.apply(format_row, axis=1)

    # Join all formatted rows with '\n\n' in between each
    SPECTRUM = '\n\n'.join(formatted_rows.tolist())

    return SPECTRUM


def CSV_convert_processing(FINAL_CSV):
    """
    Process CSV conversion using multithreading.

    :param FINAL_CSV: The list of CSV files to be converted.
    :return: The list of successfully converted CSV files.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(csv_to_msp, FINAL_CSV), total=len(FINAL_CSV), unit=" spectrums", colour="green", desc="\t  converting"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


def convert_to_msp(input_path):
    """
    Converts JSON and XML files to MSP format.

    :param input_path: The path to the directory containing the JSON and XML files.
    :return: A tuple containing the converted JSON and XML files in MSP format.
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
    xml_path = os.path.join(input_path, "XML")
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

    # CSV
    FINAL_CSV = []
    csv_to_do = False
    csv_path = os.path.join(input_path, "CSV")
    # check if there is a csv file into the directory
    for files in os.listdir(csv_path):
        if files.endswith(".csv"):
            csv_to_do = True
    if csv_to_do == True:
        time.sleep(0.02)
        print("-- CONVERTING CSV TO MSP --")
        # Concatenate all CSV to a list
        FINAL_CSV = concatenate_csv(csv_path)
        # Convert all CSV spectrum to MSP spectrum (Multithreaded)
        FINAL_CSV = CSV_convert_processing(FINAL_CSV)

    return FINAL_JSON,FINAL_XML,FINAL_CSV