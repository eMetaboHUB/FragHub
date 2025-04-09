from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import RDLogger, Chem
import pandas as pd
import globals_vars
import re
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages


def take_coresponding_mass_diff(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
                The key used in this function is "PRECURSORTYPE" which is expected to be in the dictionary.
    :return: The corresponding mass difference for the precursor type specified in the metadata dictionary.
             If the precursor type is found in the adduct_massdiff_dict, then the mass difference is returned.
             This mass difference should be a string that can be converted to a float.
             If the precursor type is not found in the adduct_massdiff_dict, None is returned.

    The adduct_massdiff_dict referenced in this function is a dictionary which has a mapping of precursor types
    to their corresponding mass differences. This dictionary should be defined globally or in the function's scope.
    """
    instrument_type = metadata_dict["INSTRUMENTTYPE"]
    if re.search("\b(IE|EI)\b", instrument_type, flags=re.IGNORECASE):
        if not metadata_dict["PRECURSORTYPE"]:
            return None

    if metadata_dict["PRECURSORTYPE"] in globals_vars.adduct_massdiff_dict:  # Check if the precursor type exists in the mass difference dictionary
        mass_diff = float(globals_vars.adduct_massdiff_dict[metadata_dict["PRECURSORTYPE"]])  # Convert the mass difference to float and assign it to a variable
        return mass_diff  # Return the mass difference
    return None  # If the precursor type was not found in the dictionary, return None

def precursor_mz_need_re_calculation(metadata_dict):
    """
    This function takes a dictionary which has metadata information as parameter.

    :param metadata_dict: A dictionary containing metadata information.
    :return: This function returns 'True' if the 'PRECURSORMZ' value in the metadata dictionary needs to be recalculated, and 'False' otherwise.
    """
    # Check if 'PRECURSORMZ' value matches the float_check_pattern
    if not re.search(globals_vars.float_check_pattern, str(metadata_dict["PRECURSORMZ"])):
        # If 'PRECURSORMZ' value doesn't match the float_check_pattern, return True indicating that the value needs to be recalculated
        return True
    elif float(re.search(globals_vars.float_check_pattern, str(metadata_dict["PRECURSORMZ"])).group(1).replace(",", ".")) <= 0.0:
        # If 'PRECURSORMZ' value is less than or equal to 0.0, return True indicating that the value needs to be recalculated.
        return True
    return False
    # If none of the above conditions are met, return False indicating that the 'PRECURSORMZ' value does not need to be recalculated.

def precursor_mz_re_calculation(spectrum, mass_diff):
    """
    Calculate the precursor m/z (mass-to-charge ratio) for a given molecular structure and a mass difference.
    :param mols: The molecular structure represented as a SMILES or InChI string.
    :param mass_diff: The mass difference to be added to the exact mass of the molecular structure.
    :return: The calculated precursor m/z value or None if the molecular structure is invalid.
    Note: This function assumes the presence of the RDKit (Chem) library for molecular structure manipulation and mass calculation.
    """
    # Initialize molecule variable as None
    mols = None

    # Assign InChi or SMILES value to the mols variable if exist in the spectrum dict
    if spectrum["INCHI"]:
        mols = spectrum["INCHI"]
    elif spectrum["SMILES"]:
        mols = spectrum["SMILES"]

    # If molecule is a string, use RDKit to turn it into a molecule object
    if isinstance(mols, str):
        mols = Chem.MolFromInchi(mols) if 'InChI=' in mols else Chem.MolFromSmiles(mols)

        # If mols is not empty, calculate the exact mass and add the mass difference
        if mols:
            exact_mass = ExactMolWt(mols)
            precursor_mz = exact_mass + mass_diff
            return precursor_mz

        # Return None if molecule string is not successfully converted to a molecule object
        return None

    # Return None if mols was not defined by any input strings in the spectrum dict
    return None

def missing_precursormz_re_calculation(metadata_dict):
    """
    :param metadata_dict: dictionary containing metadata
    :return: updated metadata dictionary
    """
    # 'PRECURSORMZ' key is checked in the metadata_dict.
    if "PRECURSORMZ" in metadata_dict:
        # The need to recalculate the 'PRECURSORMZ' is checked using the 'precursor_mz_need_re_calculation' function.
        if precursor_mz_need_re_calculation(metadata_dict):
            # If a recalculation is needed, the corresponding mass difference is taken using 'take_coresponding_mass_diff' function.
            mass_diff = take_coresponding_mass_diff(metadata_dict)
            # If mass_diff exists, the recalculation is performed and the result is updated in the metadata_dict.
            if mass_diff:
                metadata_dict["PRECURSORMZ"] = str(precursor_mz_re_calculation(metadata_dict, mass_diff))
                # If the new PRECURSORMZ value exists in the metadata_dict, the updated metadata_dict is returned.
                if metadata_dict["PRECURSORMZ"]:
                    return metadata_dict
    # If no PRECURSORMZ is present in metadata_dict or no recalculation is needed, the input metadata_dict is returned as it is.
    return metadata_dict
