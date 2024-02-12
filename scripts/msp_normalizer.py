from normalizer.values_normalizer import *
from set_parameters import parameters_dict
from peaks_filters.filters import *
import concurrent.futures
from tqdm import tqdm
import numpy as np
import re

np.set_printoptions(suppress=True)

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

def peak_list_to_np_array(peak_list, precursormz):
    """
    Converts a peak list string into a numpy array.

    :param peak_list: A string representing a peak list. Each peak is represented by a pair of values, separated by a space or colon. The first value represents the m/z (mass-to-charge ratio
    *) of the peak, and the second value represents the intensity of the peak.
    :return: A numpy array containing the peak data, with two columns for "mz" and "intensity". The "mz" column contains the m/z values, and the "intensity" column contains the corresponding
    * peak intensities.
    """
    # Convert list of tuples to numpy array
    peak_list = np.array(peak_list, dtype=float)

    # Sort the array based on the mz values
    peak_list = peak_list[peak_list[:, 0].argsort()]

    peak_list = apply_filters(peak_list, precursormz, parameters_dict)

    return peak_list

def peak_list_to_str(peak_list_np):
    """
    :param peak_list_np: The input peak list as a numpy array.
    :return: A string representation of the peak list where each row is formatted as a space-separated string of floating point values.
    """
    # Convert array to list
    peak_list_np = peak_list_np.tolist()

    # Convert list to JSON string
    peak_list_np = str(peak_list_np)

    return peak_list_np

def spectrum_cleaning(spectrum):
    """
    Parse the given spectrum to extract metadata and peak list.

    :param spectrum: The spectrum data to be parsed.
    :return: The parsed metadata with peak list.
    """
    peak_list = spectrum["PEAKS_LIST"]

    if not peak_list:
        return None

    spectrum = normalize_values(spectrum)

    if not spectrum:
        return None

    if "PRECURSORMZ" in spectrum:
        if re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])):
            spectrum["PRECURSORMZ"] = re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])).group(1)
            float_precursor_mz = float(spectrum["PRECURSORMZ"].replace(",", "."))
            if float_precursor_mz <= 0.0:
                return None
            peak_list_np = peak_list_to_np_array(peak_list, float_precursor_mz)
            if peak_list_np.size == 0:
                return None
            spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
            peak_list_np = peak_list_to_str(peak_list_np)
            spectrum["PEAKS_LIST"] = peak_list_np
            return spectrum
        else:
            return None

    return spectrum

def spectrum_cleaning_processing(spectrum_list):
    """
    :param spectrum_list: A list of spectrum data to be processed.
    :return: A list containing the results of processing the spectrum data.
    """

    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>70}".format("cleaning spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]
            results = list(executor.map(spectrum_cleaning, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    progress_bar.close()

    return final
