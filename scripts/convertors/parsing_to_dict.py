from scripts.convertors.json_to_dict import *
from scripts.convertors.msp_to_dict import *
from scripts.convertors.csv_to_dict import *
from scripts.convertors.mgf_to_dict import *
from scripts.convertors.loaders import *
import pandas as pd
import time
import json
import os
import re

def concatenate_MSP(msp_list, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
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
        spectrum_list.extend(load_spectrum_list_from_msp(files, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback))

    # Return combined list of spectra
    return spectrum_list


def concatenate_csv(csv_list, progress_callback=None, total_items_callback=None, prefix_callback=None,
                    item_type_callback=None):
    """
    Concatenates multiple CSV files into a single DataFrame, with support for progress reporting via callbacks.
    :param csv_list: List of paths to CSV files to be concatenated.
    :type csv_list: list[str]
    :param progress_callback: A function to update the progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: Concatenated DataFrame.
    :rtype: pandas.DataFrame
    """
    # Define the total number of CSV files
    total_files = len(csv_list)

    # Set total items via callback if provided
    if total_items_callback:
        total_items_callback(total_files, 0)  # total = total_files, completed = 0

    # Update the prefix dynamically via callback if provided
    if prefix_callback:
        prefix_callback("Reading CSV files:")

    # Specify the type of items being processed via callback if provided
    if item_type_callback:
        item_type_callback("csv_files")

    df_list = []  # Initialize an empty list to hold dataframes
    processed_files = 0  # Counter to track progress

    for file in csv_list:  # Iterate over all the CSV files
        # Read each CSV file as a DataFrame
        df = pd.read_csv(file, sep=";", quotechar='"', encoding="UTF-8", dtype=str)
        df['filename'] = os.path.basename(file)  # Add a 'filename' column to store the filename of each record
        df.columns = df.columns.str.lower()  # Convert all column names to lowercase
        df = df.astype(str)  # Convert all frame data to string type
        df_list.append(df)  # Append the DataFrame to the list

        # Update progress via callback if provided
        processed_files += 1
        if progress_callback:
            progress_callback(processed_files)

    # Concatenate all the DataFrames in df_list
    df = pd.concat(df_list, ignore_index=True)

    # Return the final concatenated DataFrame
    return df


def concatenate_JSON(json_list, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
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
        try:
            spectrum_list.extend(load_spectrum_list_json(files, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback))
        except:
            spectrum_list.extend(load_spectrum_list_json_2(files, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback))

    # The final list of all spectra from all JSON files is returned
    return spectrum_list

def concatenate_MGF(mgf_list, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
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
        spectrum_list.extend(load_spectrum_list_from_mgf(files, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback))
    return spectrum_list  # Return the final concatenated list of all spectra


def parsing_to_dict(input_path, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None, step_callback=None):
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
    for files in input_path:
        # Check if current file is a json file
        if files.endswith(".json"):
            # Get the full path to the json file
            json_path = files
            # Add the json file path to our list
            json_list.append(json_path)
            # Set json_to_do to True as we found at least one json file
            json_to_do = True
    # If there are json files to process
    if json_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Print a status message that the conversion has started
        if step_callback:
            step_callback("-- PARSING JSON TO DICT --")
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all json files into a single list
        FINAL_JSON = concatenate_JSON(json_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        # Reformating the JSON structure to a better structure
        FINAL_JSON = json_to_dict_processing(FINAL_JSON, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # MSP
    # Initializing an empty list to hold the final MSP data
    FINAL_MSP = []
    # Creating a list to store the paths of all MSP files found within the directory
    msp_list = []
    # A boolean variable to check if there is any MSP file to process within the directory
    msp_to_do = False
    # Looping through the directory and all subdirectories to find files
    for files in input_path:
        # Checking if the current file is an MSP file
        if files.endswith(".msp"):
            # Constructing the full path to the MSP file
            msp_path = files
            # Adding the path of the current MSP file to our list
            msp_list.append(msp_path)
            # Setting the msp_to_do flag to True as we have at least one MSP file to process
            msp_to_do = True
    # If there are any MSP files to process
    if msp_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Displaying a status message indicating the beginning of MSP to JSON conversion
        if step_callback:
            step_callback("-- PARSING MSP TO DICT --")
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all found MSP files into a single list
        FINAL_MSP = concatenate_MSP(msp_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        # Converting each MSP spectrum to a JSON spectrum
        FINAL_MSP = msp_to_dict_processing(FINAL_MSP, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # MGF
    # Initializing an empty list to hold the final MGF (Mascot Generic Format) data
    FINAL_MGF = []
    # A list to store the paths of all MGF files found within the directory
    mgf_list = []
    # A boolean variable to check if there is any MGF file to process within the directory
    mgf_to_do = False
    # Looping over each directory, subdirectory, and file in the provided directory
    for files in input_path:
        # Checking if the file has a .mgf extension
        if files.endswith(".mgf"):
            # Constructing the full path to the MGF file
            mgf_path = files
            # Adding the path of the current MGF file to our list
            mgf_list.append(mgf_path)
            # Setting the mgf_to_do flag to True as we have at least one MGF file to process
            mgf_to_do = True
    # Only proceed if there are MGF files to process
    if mgf_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Printing a status message signaling the start of the MGF to JSON conversion
        if step_callback:
            step_callback("-- PARSING MGF TO DICT --")
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenating all found MGF files into a single list
        FINAL_MGF = concatenate_MGF(mgf_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        # Converting each MGF spectrum to a JSON spectrum
        FINAL_MGF = mgf_to_dict_processing(FINAL_MGF, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # CSV
    # Initializing an empty list to contain the final CSV data
    FINAL_CSV = []
    # A list to keep track of all CSV file paths found within the directory
    csv_list = []
    # A boolean flag to check if there is any CSV file to process
    csv_to_do = False
    # Walking through the directory and all its subdirectories to look for files
    for files in input_path:
        # Checking if the current file has the extension '.csv'
        if files.endswith(".csv"):
            # If it's a CSV file, getting full path to it
            csv_path = files
            # Adding this CSV file path to our list
            csv_list.append(csv_path)
            # Setting the csv_to_do flag to True as there's at least one CSV file to process
            csv_to_do = True
    # Go ahead only if there are CSV files to process
    if csv_to_do == True:
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Printing a status message to signify the start of CSV to JSON conversion
        if step_callback:
            step_callback("-- PARSING CSV TO DICT --")
        # Sleep for a short time to correctly display progress bar
        time.sleep(0.01)
        # Concatenate all CSV files into one list
        FINAL_CSV = concatenate_csv(csv_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        # Convert the CSV data to JSON
        FINAL_CSV = csv_to_dict_processing(FINAL_CSV, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    return FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF
