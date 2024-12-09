from .json_to_dict import *
from .msp_to_dict import *
from .xml_to_json import *
from .csv_to_dict import *
from .mgf_to_dict import *
from .loaders import *
from tqdm import tqdm
import pandas as pd
import time
import json
import os
import re

def concatenate_MSP(msp_list):
    """
    This function concatenates multiple MSP files into a single spectrum list.

    :param msp_list: A list of paths to MSP files. List structure: [file1, file2, file3, ...]
    :return: Returns a list of spectra loaded from the MSP files.

    The function loads the MSP files in the order they appear in the msp_list, combining
    them into one continuous list of spectra.
    """
    # Initialize empty list to store spectra
    spectrum_list = []

    # Loop over each file in the msp_list
    for files in msp_list:
        # The extend method appends elements from the iterable (the result of load_spectrum_list_from_msp function)
        # to the end of the list.
        spectrum_list.extend(load_spectrum_list_from_msp(files))

    # Return combined list of spectra
    return spectrum_list

def concatenate_xml(xml_list):
    """
    Concatenates the XML contents of a list of XML files and adds a filename tag to each XML.

    :param xml_list: A list of XML file paths.
    :return: A list containing the final XML contents.
    """
    FINAL_XML = []

    pbar = tqdm(xml_list, total=len(xml_list), unit=" spectrums", colour="green", desc="{:>70}".format("concatenate"))
    for files in pbar:
        if files.endswith(".xml"):
            file_name = os.path.basename(files.replace(".xml", ""))
            with open(files, "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            xml_content = re.sub("</id>\n", f"</id>\n  <filename>{file_name}</filename>\n", xml_content)
            FINAL_XML.extend([xml_content])
    pbar.close()
    return FINAL_XML

def concatenate_csv(csv_list):
    """
    Concatenates multiple CSV files into a single DataFrame.
    :param csv_list: List of paths to CSV files to be concatenated.
    :type csv_list: list[str]
    :return: Concatenated DataFrame.
    :rtype: pandas.DataFrame
    """
    df_list = []  # Initialize an empty list to hold dataframes
    for file in csv_list:  # Iterate over all the csv files
        # Read each csv file as a dataframe
        df = pd.read_csv(file, sep=";", quotechar='"', encoding="UTF-8", dtype=str)
        df['filename'] = os.path.basename(file)  # Add a 'filename' column to store the filename of each record
        df.columns = df.columns.str.lower()  # Convert all column names to lowercase
        df = df.astype(str)  # Convert all frame data to string data type
        df_list.append(df)  # Append the dataframe to the list
    # Concatenate all the dataframes in df_list
    df = pd.concat(df_list, ignore_index=True)
    return df  # Return the final concatenated dataframe

def concatenate_JSON(json_list):
    """
    This function is used to concatenate JSON files into a single list.

    Parameters:
        json_list (list[str]) : A list of JSON file paths.

    Returns:
        list[dict] : The concatenated list with spectrum dictionaries.

    Usage:
        The function is used when there are multiple JSON files, but the application
        needs them in one collective list. Used mainly in file conversion or
        preprocessing before data analysis.
    """

    # Initialize an empty list to keep all spectrums
    spectrum_list = []

    # Loop over each file in the json_list
    for files in json_list:
        # Using extend method to add elements of individual spectrums
        # from each file to the global spectrum_list
        # 'load_spectrum_list_json()' is used to load list of spectrum
        # from individual json file
        spectrum_list.extend(load_spectrum_list_json(files))

    # The final list of all spectra from all JSON files is returned
    return spectrum_list

def concatenate_MGF(mgf_list):
    """
    Concatenates multiple MGF files into a single spectrum list.

    This function is instrumental in the conversion of MGF files as it binds
    all MGF files' data into a single list by iterating over the provided list of file paths.
    Once each file's spectrum list has been loaded via the external function `load_spectrum_list_from_mgf()`,
    it extends the main list `spectrum_list` with the new spectrum data from each file.
    This way, all the spectral data from all MGF files are consolidated into one list.

    :param mgf_list: A list of paths to MGF files. These are the files that the function will concatenate
                     the spectral data from. Each item in the list should be a string representing an absolute or
                     relative path to a MGF file.
    :type mgf_list: list[str])
    :return: The concatenated spectrum list. This list contains all the spectra from all MGF files provided.
    :rtype: list[str]
    """
    spectrum_list = []  # Initialize an empty list to store spectra
    for files in mgf_list:  # Loop over each file path in the mgf_list
        # Extend the spectrum_list with spectral data from current file
        # The load_spectrum_list_from_mgf() is supposed to return a list of spectra from a single MGF file.
        spectrum_list.extend(load_spectrum_list_from_mgf(files))
    return spectrum_list  # Return the final concatenated list of all spectra

def parsing_to_dict(input_path):
    """
    This function converts file data (JSON, XML, CSV, MSP, MGF) in a directory to JSON format.
    :param input_path: The path of the input directory where the files are located.
    :return: A tuple containing the converted data from each file type to JSON format. The order of the elements in the tuple is as follows:
        - FINAL_MSP: The converted data from MSP files to JSON format.
        - FINAL_XML: The converted data from XML files to JSON format.
        - FINAL_CSV: The converted data from CSV files to JSON format.
        - FINAL_JSON: The converted data from JSON files to JSON format (not actually converted, just collected).
        - FINAL_MGF: The converted data from MGF files to JSON format.
    """
    # JSON
    # Initializing an empty list to contain the final JSON data
    FINAL_JSON = []
    # A list to store all the paths of json files found in the directory
    json_list = []
    # A boolean variable to check if there is any json file in the directory
    json_to_do = False
    # Iterate over all directories, subdirectories and files in the given directory
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Check if current file is a json file
            if file.endswith(".json"):
                # Get the full path to the json file
                json_path = os.path.join(root, file)
                # Add the json file path to our list
                json_list.append(json_path)
                # Set json_to_do to True as we found at least one json file
                json_to_do = True
    # If there are json files to process
    if json_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Print a status message that the conversion has started
        print("{:>70}".format("-- PARSING JSON TO DICT --"))
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all json files into a single list
        FINAL_JSON = concatenate_JSON(json_list)
        # Reformating the JSON structure to a better structure
        FINAL_JSON = json_to_dict_processing(FINAL_JSON)

    # MSP
    # Initializing an empty list to hold the final MSP data
    FINAL_MSP = []
    # Creating a list to store the paths of all MSP files found within the directory
    msp_list = []
    # A boolean variable to check if there is any MSP file to process within the directory
    msp_to_do = False
    # Looping through the directory and all subdirectories to find files
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Checking if the current file is an MSP file
            if file.endswith(".msp"):
                # Constructing the full path to the MSP file
                msp_path = os.path.join(root, file)
                # Adding the path of the current MSP file to our list
                msp_list.append(msp_path)
                # Setting the msp_to_do flag to True as we have at least one MSP file to process
                msp_to_do = True
    # If there are any MSP files to process
    if msp_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Displaying a status message indicating the beginning of MSP to JSON conversion
        print("{:>70}".format("-- PARSING MSP TO DICT --"))
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all found MSP files into a single list
        FINAL_MSP = concatenate_MSP(msp_list)
        # Converting each MSP spectrum to a JSON spectrum
        FINAL_MSP = msp_to_dict_processing(FINAL_MSP)

    # MGF
    # Initializing an empty list to hold the final MGF (Mascot Generic Format) data
    FINAL_MGF = []
    # A list to store the paths of all MGF files found within the directory
    mgf_list = []
    # A boolean variable to check if there is any MGF file to process within the directory
    mgf_to_do = False
    # Looping over each directory, subdirectory, and file in the provided directory
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Checking if the file has a .mgf extension
            if file.endswith(".mgf"):
                # Constructing the full path to the MGF file
                mgf_path = os.path.join(root, file)
                # Adding the path of the current MGF file to our list
                mgf_list.append(mgf_path)
                # Setting the mgf_to_do flag to True as we have at least one MGF file to process
                mgf_to_do = True
    # Only proceed if there are MGF files to process
    if mgf_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Printing a status message signaling the start of the MGF to JSON conversion
        print("{:>70}".format("-- PARSING MGF TO DICT --"))
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all found MGF files into a single list
        FINAL_MGF = concatenate_MGF(mgf_list)
        # Converting each MGF spectrum to a JSON spectrum
        FINAL_MGF = mgf_to_dict_processing(FINAL_MGF)

    # XML
    # Initialize an empty list to contain the final XML data
    FINAL_XML = []
    # A list to store the paths of all xml files found in the directory
    xml_list = []
    # A boolean variable to check if there is any xml file in the directory
    xml_to_do = False
    # Loop over all directories, subdirectories, and files in the provided directory
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Check if the current file is an xml file
            if file.endswith(".xml"):
                # Get the full path of the xml file
                xml_path = os.path.join(root, file)
                # Append the xml file path to the list
                xml_list.append(xml_path)
                # Set xml_to_do to True since we have at least one xml file
                xml_to_do = True
    # If there are xml files to be processed
    if xml_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Print a status message that the conversion process has begun
        print("{:>70}".format("-- CONVERTING XML TO JSON --"))
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all xml files into one list
        FINAL_XML = concatenate_xml(xml_list)
        # Formatting the XML structure into JSON structure
        FINAL_XML = xml_to_json_processing(FINAL_XML)

    # CSV
    # Initializing an empty list to contain the final CSV data
    FINAL_CSV = []
    # A list to keep track of all CSV file paths found within the directory
    csv_list = []
    # A boolean flag to check if there is any CSV file to process
    csv_to_do = False
    # Walking through the directory and all its subdirectories to look for files
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Checking if the current file has the extension '.csv'
            if file.endswith(".csv"):
                # If it's a CSV file, getting full path to it
                csv_path = os.path.join(root, file)
                # Adding this CSV file path to our list
                csv_list.append(csv_path)
                # Setting the csv_to_do flag to True as there's at least one CSV file to process
                csv_to_do = True
    # Go ahead only if there are CSV files to process
    if csv_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Printing a status message to signify the start of CSV to JSON conversion
        print("{:>70}".format("-- PARSING CSV TO DICT --"))
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenate all CSV files into one list
        FINAL_CSV = concatenate_csv(csv_list)
        # Convert the CSV data to JSON
        FINAL_CSV = csv_to_dict_processing(FINAL_CSV)

    return FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF
