from rdkit.Chem.Descriptors import ExactMolWt
from values_normalizer import *
from ast import literal_eval
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

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global adduct_massdiff_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_massdiff_dict = dict(zip(adduct_dataframe['fraghub_default'], adduct_dataframe['massdiff']))

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

    peak_list = apply_filters(peak_list, precursormz)

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

def precursor_mz_re_calculation(spectrum, mass_diff):
    """
    Calculate the precursor m/z (mass-to-charge ratio) for a given molecular structure and a mass difference.

    :param mols: The molecular structure represented as a SMILES or InChI string.
    :param mass_diff: The mass difference to be added to the exact mass of the molecular structure.
    :return: The calculated precursor m/z value or None if the molecular structure is invalid.

    Note: This function assumes the presence of the RDKit (Chem) library for molecular structure manipulation and mass calculation.
    """
    mols = None
    if spectrum["INCHI"]:
        mols = spectrum["INCHI"]
    elif spectrum["SMILES"]:
        mols = spectrum["SMILES"]

    if isinstance(mols, str):
        mols = Chem.MolFromInchi(mols) if 'InChI=' in mols else Chem.MolFromSmiles(mols)
        if mols:
            exact_mass = ExactMolWt(mols)
            precursor_mz = exact_mass + mass_diff

            return precursor_mz
        return None

    return None

def spectrum_cleaning(spectrum):
    """
    :param spectrum: a dictionary containing information about a spectrum
    :return: the cleaned spectrum dictionary or None if cleaning could not be performed

    This method performs cleaning on a spectrum by normalizing values, recalculating the precursor mass if necessary,
    and converting the peak list format.

    The `spectrum` parameter should be a dictionary with the following keys:
    - "PEAKS_LIST": a list of peaks
    - "PRECURSORMZ": the precursor mass-to-charge ratio
    - "PRECURSORTYPE": the type of precursor ion
    - "INCHI": the InChI representation of the molecule
    - "SMILES": the SMILES representation of the molecule

    If the `PEAKS_LIST` key does not exist or is empty, the method returns None.

    If normalization fails or the normalized spectrum is empty, the method returns None.

    If the "PRECURSORMZ" key exists and the value is a valid floating-point number, the method performs the following steps:
    1. Extracts the float value from the string representation of the precursor mass.
    2. Checks if the float value is less than or equal to 0.0. If true, returns None.
    3. Converts the peak list to a numpy array using the float value as the precursor mass-to-charge ratio.
    4. Checks if the peak list numpy array is empty. If true, returns None.
    5. Updates the "NUM PEAKS" key in the spectrum dictionary to the number of peaks in the peak list numpy array.
    6. Converts the peak list numpy array to a string representation and updates the "PEAKS_LIST" key in the spectrum dictionary.

    If the "PRECURSORMZ" key does not exist or the value is not a valid floating-point number, the method performs the following steps:
    1. Checks if the "PRECURSORTYPE" key exists and has a corresponding mass difference value in the `adduct_massdiff_dict` dictionary.
    2. If true, retrieves the mass difference for the precursor type.
    3. Checks if the "INCHI" or "SMILES" key exists in the spectrum dictionary.
    4. If true, recalculates the precursor mass using the molecule representation and the mass difference.
    5. Checks if the recalculated precursor mass is a valid floating-point number. If not, returns None.
    6. Converts the peak list to a numpy array using the recalculated precursor mass.
    7. Checks if the peak list numpy array is empty. If true, returns None.
    8. Updates the "NUM PEAKS" key in the spectrum dictionary to the number of peaks in the peak list numpy array.
    9. Converts the peak list numpy array to a string representation and updates the "PEAKS_LIST" key in the spectrum dictionary.

    If none of the above conditions are met, the method returns None.

    Example usage:
    cleaned_spectrum = spectrum_cleaning(spectrum)
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
                if "PRECURSORTYPE" in spectrum:
                    if spectrum["PRECURSORTYPE"] in adduct_massdiff_dict:
                        mass_diff = float(adduct_massdiff_dict[spectrum["PRECURSORTYPE"]])
                        spectrum["PRECURSORMZ"] = precursor_mz_re_calculation(spectrum, mass_diff)
                        if spectrum["PRECURSORMZ"]:
                            if re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])):
                                spectrum["PRECURSORMZ"] = re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])).group(1)
                                float_precursor_mz = float(spectrum["PRECURSORMZ"].replace(",", "."))
                                if float_precursor_mz <= 0.0:  # normalement impossible dans ce else puisqu'il a été recalculer
                                    return None
                                peak_list_np = peak_list_to_np_array(peak_list, float_precursor_mz)
                                if peak_list_np.size == 0:
                                    return None
                                spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
                                peak_list_np = peak_list_to_str(peak_list_np)
                                spectrum["PEAKS_LIST"] = peak_list_np
                                return spectrum
                return None
            peak_list_np = peak_list_to_np_array(peak_list, float_precursor_mz)
            if peak_list_np.size == 0:
                return None
            spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
            peak_list_np = peak_list_to_str(peak_list_np)
            spectrum["PEAKS_LIST"] = peak_list_np
            return spectrum
        else:
            # re calculer la masse de l'ion precurseur
            if "PRECURSORTYPE" in spectrum:
                if spectrum["PRECURSORTYPE"] in adduct_massdiff_dict:
                    mass_diff = float(adduct_massdiff_dict[spectrum["PRECURSORTYPE"]])
                    spectrum["PRECURSORMZ"] = precursor_mz_re_calculation(spectrum, mass_diff)
                    if spectrum["PRECURSORMZ"]:
                        if re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])):
                            spectrum["PRECURSORMZ"] = re.search(float_check_pattern, str(spectrum["PRECURSORMZ"])).group(1)
                            float_precursor_mz = float(spectrum["PRECURSORMZ"].replace(",", "."))
                            if float_precursor_mz <= 0.0: # normalement impossible dans ce else puisqu'il a été recalculer
                                return None
                            peak_list_np = peak_list_to_np_array(peak_list, float_precursor_mz)
                            if peak_list_np.size == 0:
                                return None
                            spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
                            peak_list_np = peak_list_to_str(peak_list_np)
                            spectrum["PEAKS_LIST"] = peak_list_np
                            return spectrum

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
