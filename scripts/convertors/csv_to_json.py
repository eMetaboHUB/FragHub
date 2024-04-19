from .keys_convertor import *
from tqdm import tqdm
import pandas as pd
import re

global peak_list_csv_to_json_pattern
peak_list_csv_to_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def parse_peak_list(peak_list_string):
    """
    Parses a peak list string into a list of peak values.
    :param peak_list_string: A string representing the peak list.
    :return: A list of peak values, where each peak value is a list containing the m/z value and intensity.
    """
    # Here, 're.findall' method is being used to find all the occurrences of peak_list_csv_to_json_pattern in peak_list_string.
    peaks = re.findall(peak_list_csv_to_json_pattern, peak_list_string)

    # Peaks are being converted to a list of two values' lists, where each sub-list contains the float values of m/z and intensity.
    return [[float(mz), float(intensity)] for mz, intensity in peaks]

def csv_to_json_processing(FINAL_CSV):
    """
    Convert a CSV file to a list of JSON objects.

    :param FINAL_CSV: A pandas DataFrame representing the CSV file.
    :return: A list of JSON objects converted from the CSV file.
    """
    # Here, 'FINAL_CSV.to_dict()' method is being used to convert each row of CSV into a dictionary and all dictionaries are combined into a list.
    json_list = FINAL_CSV.to_dict('records')

    # The 'if' and 'elif' checks whether 'peaks' or 'peaks_list' is in the row, and if so, calls the parse_peak_list function for them.
    for row in json_list:
        if "peaks" in row:
            row["peaks"] = parse_peak_list(row["peaks"])
        elif "peaks_list" in row:
            row["peaks"] = parse_peak_list(row["peaks_list"])

    # A progress bar is created to provide visual output to the user about the conversion process progression.
    with tqdm(total=len(json_list), unit=" spectrums", colour="green", desc="{:>70}".format("converting CSV spectrums")) as pbar:
        for i in range(len(json_list)):
            # Keys of each JSON object are being converted by calling 'convert_keys' function for each JSON object.
            json_list[i] = convert_keys(json_list[i])

            # The progress bar is being updated after converting each JSON object.
            pbar.update()
        pbar.close()

    # Finally, the converted list of JSON objects is being returned.
    return json_list

