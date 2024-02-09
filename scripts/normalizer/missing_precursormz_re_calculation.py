from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import RDLogger, Chem
import pandas as pd
import re
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global adduct_massdiff_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_massdiff_dict = dict(zip(adduct_dataframe['fraghub_default'], adduct_dataframe['massdiff']))

def take_coresponding_mass_diff(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The corresponding mass difference for the precursor type specified in the metadata dictionary. If the precursor type is found in the adduct_massdiff_dict, the mass difference
    * is returned. Otherwise, None is returned.
    """
    if metadata_dict["PRECURSORTYPE"] in adduct_massdiff_dict:
        mass_diff = float(adduct_massdiff_dict[metadata_dict["PRECURSORTYPE"]])
        return mass_diff

    return None

def precursor_mz_need_re_calculation(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: True if the precursor m/z value in the metadata dictionary needs to be re-calculated, False otherwise.

    This method checks whether the precursor m/z value in the given metadata dictionary needs to be re-calculated. If the precursor m/z value does not match the float_check_pattern or if
    * it is less than or equal to 0.0, the method returns True indicating that the value needs to be re-calculated. Otherwise, it returns False.
    """
    if not re.search(float_check_pattern, str(metadata_dict["PRECURSORMZ"])):
        return True
    elif float(re.search(float_check_pattern, str(metadata_dict["PRECURSORMZ"])).group(1).replace(",", ".")) <= 0.0:
        return True

    return False

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

def missing_precursormz_re_calculation(metadata_dict):
    """
    :param metadata_dict: dictionary containing metadata
    :return: updated metadata dictionary
    """
    if "PRECURSORMZ" in metadata_dict:
        if precursor_mz_need_re_calculation(metadata_dict):
            mass_diff = take_coresponding_mass_diff(metadata_dict)
            if mass_diff:
                metadata_dict["PRECURSORMZ"] = str(precursor_mz_re_calculation(metadata_dict, mass_diff))
                if metadata_dict["PRECURSORMZ"]:
                    return metadata_dict

    return metadata_dict