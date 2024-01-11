from values_normalizer import *
import concurrent.futures
from tqdm import tqdm
from filters import *
import pandas as pd
import numpy as np
import ijson
import json
import os
import re

np.set_printoptions(suppress=True)

global keys_dict
Key_dataframe = pd.read_csv(os.path.abspath("../datas/key_to_convert.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
keys_dict = dict(zip(Key_dataframe['known_synonym'], Key_dataframe['fraghub_default'].str.upper()))

global keys_list
keys_list = ['FILENAME',
             'PREDICTED',
             'FRAGHUBID',
             'SPECTRUMID',
             'RESOLUTION',
             'SYNON',
             'CHARGE',
             'IONIZATION',
             'MSLEVEL',
             'FRAGMENTATIONMODE',
             'NAME',
             'PRECURSORMZ',
             'EXACTMASS',
             'AVERAGEMASS',
             'PRECURSORTYPE',
             'INSTRUMENTTYPE',
             'INSTRUMENT',
             'SMILES',
             'INCHI',
             'INCHIKEY',
             'COLLISIONENERGY',
             'FORMULA',
             'RETENTIONTIME',
             'IONMODE',
             'COMMENT',
             'NUM PEAKS',
             'PEAKS_LIST']

def convert_keys(metadata_dict):
    """
    Convert keys in metadata_dict based on the provided keys_dict and keys_list.

    :param metadata_dict: A dictionary containing metadata information.
    :return: A dictionary with converted keys based on the provided keys_dict and keys_list.
    """

    converted = {keys_dict[key]: val for key, val in metadata_dict.items() if
                 key in keys_dict and keys_dict[key] in keys_list}

    converted.update({key: "" for key in keys_list if key not in converted})

    return converted

def peak_list_to_np_array(peak_list, precursormz):
    """
    Converts a peak list string into a numpy array.

    :param peak_list: A string representing a peak list. Each peak is represented by a pair of values, separated by a space or colon. The first value represents the m/z (mass-to-charge ratio
    *) of the peak, and the second value represents the intensity of the peak.
    :return: A numpy array containing the peak data, with two columns for "mz" and "intensity". The "mz" column contains the m/z values, and the "intensity" column contains the corresponding
    * peak intensities.
    """
    peak_list = json.loads(peak_list)

    # Convert list of tuples to numpy array
    peak_list = np.array(peak_list, dtype=float)

    # Sort the array based on the mz values
    peak_list = peak_list[peak_list[:, 0].argsort()]

    peak_list = apply_filters(peak_list, precursormz)

    if peak_list.size == 0:
        return ''

    return peak_list

def peak_list_to_str(peak_list_np):
    """
    :param peak_list_np: The input peak list as a numpy array.
    :return: A string representation of the peak list where each row is formatted as a space-separated string of floating point values.
    """
    # Convertir l'array en liste
    peak_list_np = peak_list_np.tolist()

    # Convertir la liste en chaîne JSON
    peak_list_np = json.dumps(peak_list_np)

    return peak_list_np

def msp_parser(spectrum):
    """
    Parse the given spectrum to extract metadata and peak list.

    :param spectrum: The spectrum data to be parsed.
    :return: The parsed metadata with peak list.
    """
    peaks_list = spectrum["PEAKS_LIST"]

    spectrum = convert_keys(spectrum)
    spectrum = normalize_values(spectrum)

    if peaks_list is not None:
        spectrum["PEAKS_LIST"] = peaks_list

    if not spectrum:
        return {}

    if "PRECURSORMZ" in spectrum:
        if spectrum["PRECURSORMZ"]:
            try:
                peak_list_np = peak_list_to_np_array(spectrum["PEAKS_LIST"], float(spectrum["PRECURSORMZ"].replace(",", ".")))
                if not peak_list_np:
                    return {}
                peak_list_np = peak_list_to_str(peak_list_np)
                spectrum["PEAKS_LIST"] = peak_list_np
                return spectrum
            except:
                return {}

    return spectrum

def msp_cleaning_processing(spectrum_list):
    """
    :param spectrum_list: A list of spectrum data to be processed.
    :return: A list containing the results of processing the spectrum data.
    """

    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>80}".format("cleaning file"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]
            results = list(executor.map(msp_parser, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final
