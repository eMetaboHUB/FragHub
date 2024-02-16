from tqdm import tqdm
from .keys_convertor import *
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
    peaks = re.findall(peak_list_csv_to_json_pattern, peak_list_string)

    return [[float(mz), float(intensity)] for mz, intensity in peaks]

def csv_to_json_processing(FINAL_CSV):
    """
    Convert a CSV file to a list of JSON objects.

    :param FINAL_CSV: A pandas DataFrame representing the CSV file.
    :return: A list of JSON objects converted from the CSV file.
    """
    json_list = FINAL_CSV.to_dict('records')

    for row in json_list:
        if "peaks" in row:
            row["peaks"] = parse_peak_list(row["peaks"])
        elif "PEAKS_LIST" in row:
            row["PEAKS_LIST"] = parse_peak_list(row["PEAKS_LIST"])

    with tqdm(total=len(json_list), unit=" spectrums", colour="green", desc="{:>70}".format("converting CSV spectrums")) as pbar:
        for i in range(len(json_list)):
            json_list[i] = convert_keys(json_list[i])
            pbar.update()

    return json_list
