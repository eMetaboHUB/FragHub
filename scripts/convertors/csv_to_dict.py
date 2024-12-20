from .keys_convertor import *
import pandas as pd
import globals_vars
import re

def parse_peak_list(peak_list_string):
    """
    Parses a peak list string into a list of peak values.
    :param peak_list_string: A string representing the peak list.
    :return: A list of peak values, where each peak value is a list containing the m/z value and intensity.
    """
    # Here, 're.findall' method is being used to find all the occurrences of peak_list_csv_to_json_pattern in peak_list_string.
    peaks = re.findall(globals_vars.peak_list_json_pattern, peak_list_string)
    del peak_list_string

    # Peaks are being converted to a list of two values' lists, where each sub-list contains the float values of m/z and intensity.
    return [[float(mz), float(intensity)] for mz, intensity in peaks]

def csv_to_dict_processing(FINAL_CSV, progress_callback=None, total_items_callback=None, prefix_callback=None,
                           item_type_callback=None):
    """
    Convert a CSV file (as a DataFrame) to a list of JSON objects, with progress reporting via callbacks.

    :param FINAL_CSV: A pandas DataFrame representing the CSV file.
    :param progress_callback: A function to update progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: A list of JSON objects converted from the CSV file.
    """
    # Convert the DataFrame to a list of dictionaries (JSON-like objects)
    FINAL_CSV = FINAL_CSV.to_dict('records')  # Each row of CSV is converted into a dictionary

    # Set total items via callback if provided
    if total_items_callback:
        total_items_callback(len(FINAL_CSV), 0)  # total = len(FINAL_CSV), completed = 0

    # Update the prefix dynamically via callback if provided
    if prefix_callback:
        prefix_callback("Parsing CSV spectrums:")

    # Specify the type of items being processed via callback if provided
    if item_type_callback:
        item_type_callback("spectra")

    # Process progress
    processed_items = 0

    # Iterate over each row to process 'peaks' or 'peaks_list' and convert keys
    for i in range(len(FINAL_CSV)):
        row = FINAL_CSV[i]

        # Process `peaks` and `peaks_list` if they exist in the row
        if "peaks" in row:
            row["peaks"] = parse_peak_list(row["peaks"])
        elif "peaks_list" in row:
            row["peaks"] = parse_peak_list(row["peaks_list"])

        # Convert keys of the row
        FINAL_CSV[i] = convert_keys(row)

        # Update progress via callback if provided
        processed_items += 1
        if progress_callback:
            progress_callback(processed_items)

    # The final processed list is returned
    return FINAL_CSV


