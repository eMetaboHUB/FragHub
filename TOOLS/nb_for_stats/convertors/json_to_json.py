from .keys_convertor import *
import concurrent.futures
from tqdm import tqdm
import re

global peak_list_json_to_json_pattern
peak_list_json_to_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")


def parse_MoNA_peak_list(peak_list_string):
    """
    Parse a peak list string from the MoNA database into a list of peaks.

    :param peak_list_string: A string containing the peak list data in JSON format.
    :return: A list of peaks, where each peak is represented as a list containing the m/z value and intensity.
    :rtype: list[list[float, float]]
    """
    peaks = re.findall(peak_list_json_to_json_pattern, peak_list_string)

    return [[float(mz), float(intensity)] for mz, intensity in peaks]

def convert_MoNA_json(json_dict):
    """
    :convert_MoNA_json:

    Converts a JSON dictionary representing a MoNA (MassBank of North America) entry into a formatted dictionary.

    :param json_dict: A dictionary representing a MoNA entry in JSON format.
    :return: A dictionary containing the converted information from the MoNA entry.

    """
    dict_final = {}

    try:
        dict_final["compound_name"] = json_dict["compound"][0]["names"][0]["name"]
    except:
        dict_final["compound_name"] = ''

    try:
        for i in range(len(json_dict["compound"][0]["metaData"])):
            if json_dict["compound"][0]["metaData"][i]["name"] in ["molecular formula", "SMILES", "InChI", "InChIKey"] and not json_dict["compound"][0]["metaData"][i]["computed"]:
                dict_final[json_dict["compound"][0]["metaData"][i]["name"].lower()] = json_dict["compound"][0]["metaData"][i]["value"]
    except:
        pass

    try:
        dict_final["spectrum_id"] = json_dict["id"]
    except:
        pass

    try:
        for i in range(len(json_dict["metaData"])):
            if json_dict["metaData"][i]["name"] in ["instrument", "instrument type", "ms level", "ionization", "retention time", "ionization mode", "precursor type", "collision energy", "precursor m/z"]:
                dict_final[json_dict["metaData"][i]["name"].lower()] = json_dict["metaData"][i]["value"]
    except:
        pass

    try:
        dict_final["filename"] = json_dict["filename"]
    except:
        pass

    try:
        peak_list_string = json_dict["spectrum"]
        dict_final["peaks"] = parse_MoNA_peak_list(peak_list_string)
    except:
        pass

    return dict_final

def parse_others_json_peak_list(peak_list):
    """
    Parses a JSON peak list and returns a list of peaks.

    :param peak_list: A JSON string representing the peak list.
    :return: A list of peak tuples where each tuple contains the mz (float) and intensity (float) values.
    """
    peak_list = re.findall(peak_list_json_to_json_pattern, peak_list)

    return [[float(mz), float(intensity)] for mz, intensity in peak_list]

def json_to_json(json_dict):
    """
    :param json_dict: A dictionary representing a JSON object.
    :return: The converted JSON object.

    This method takes a JSON object represented as a dictionary and checks if it contains the keys "compound", "id", "metaData", and "spectrum". If all of these keys are present, the method
    * calls the function "convert_MoNA_json" to convert the JSON object. Otherwise, it returns the original JSON object.

    Note: The "convert_MoNA_json" function is not implemented in this code, so you need to define it separately before using this method.
    """

    keys_to_check = ["compound", "id", "metaData", "spectrum", "filename"]
    peak_list_keys = ["spectrum", "peaks_json", "peaks"]

    if all(key in json_dict for key in keys_to_check):
        json_dict = convert_MoNA_json(json_dict)
        json_dict = convert_keys(json_dict)
        return json_dict
    else:
        for key in peak_list_keys:
            if key in json_dict:
                json_dict[key] = parse_others_json_peak_list(str(json_dict[key]))
                json_dict = convert_keys(json_dict)
                return json_dict
    return None

def json_to_json_processing(FINAL_JSON):
    """
    Converts a list of JSON objects into processed JSON data.

    :param FINAL_JSON: A list of JSON objects to be processed.
    :return: A list of processed JSON data.

    This method takes a list of JSON objects and processes them using the 'json_to_json' method. The list is divided into chunks of size 'chunk_size'. Each chunk is processed concurrently
    * using a ThreadPoolExecutor, and the results are collected and returned as a list of processed JSON data.

    Example usage:
    ```
    input_json = [...] # List of JSON objects
    processed_json = json_to_json_processing(input_json)
    ```
    """
    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(FINAL_JSON), unit=" spectrums", colour="green", desc="{:>70}".format("converting JSON spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(FINAL_JSON), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = FINAL_JSON[i:i + chunk_size]
            results = list(executor.map(json_to_json, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final