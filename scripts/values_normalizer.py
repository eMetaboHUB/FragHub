from fuzzywuzzy import process
from mols_calculation import *
import numpy as np
import itertools
import re
import os

global repair_inchi_pattern
repair_inchi_pattern = re.compile(r"^(inchi=)?", flags=re.IGNORECASE)

global smiles_pattern
smiles_pattern = re.compile(r"[^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,}", flags=re.IGNORECASE) # Match smiles

global inchi_pattern
inchi_pattern = re.compile(r"InChI=.*|\/[0-9A-Z]*\/", flags=re.IGNORECASE) # Match inchi

global inchikey_pattern
inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey

global adduct_pattern
adduct_pattern = re.compile(r"(\[([A-Za-z0-9\+\-\(\)]*)\]((?:[0-9]*)?[\+\-\*])*)(?:\/|$)?")

global charge_pattern
charge_pattern = re.compile(r"((\+|\-|^)\d+)|(\+|\-)")

global adduct_pattern_2
adduct_pattern_2 = re.compile(r"([A-Za-z0-9\+\-\(\)\*]*)")

global ending_by_charge_pattern
ending_by_charge_pattern = re.compile(r"(\d)?([\+\-\*])$")

global ionmode_pos_pattern
ionmode_pos_pattern = re.compile(r"^p|^\+|^pos", flags=re.IGNORECASE)

global ionmode_neg_pattern
ionmode_neg_pattern = re.compile(r"^n|^\-|^neg", flags=re.IGNORECASE)

global precursortype_pos_pattern
precursortype_pos_pattern = re.compile(r"\][\+\*]*$")

global precursortype_neg_pattern
precursortype_neg_pattern = re.compile(r"\][\-\*]*$")

global ms_level_pattern
ms_level_pattern = re.compile(r"(?:ms)?(\d)", flags=re.IGNORECASE)

global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)

global empty_pattern
empty_pattern = re.compile(r"(^CCS:( .*)?)|(^\$:00in-source( .*)?)|(^0( .*)?)|(^0\.0( .*)?)|(^$)|(^na( .*)?)|(^n/a( .*)?)|(^nan( .*)?)|(^unknown( .*)?)|(^unknow( .*)?)|(^none( .*)?)|(^\?( .*)?)|(^unk( .*)?)|(^x( .*)?)", flags=re.IGNORECASE)

global adduct_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_dict = dict(zip(adduct_dataframe['known_adduct'], adduct_dataframe['fraghub_default']))

global adduct_massdiff_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_massdiff_dict = dict(zip(adduct_dataframe['fraghub_default'], adduct_dataframe['massdiff']))

global instruments_list
instruments_list = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instrument_types_list = instruments_list['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_list = instruments_list['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_list = [''.join(pair) for pair in itertools.zip_longest(instruments_list, instrument_types_list, fillvalue='')]

global instruments_dict
instruments_dict = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instruments_dict['INIT_INSTRUMENT'] = instruments_dict['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_dict['INIT_INSTRUMENT_TYPE'] = instruments_dict['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_dict['INDEX'] = instruments_dict['INIT_INSTRUMENT'] + instruments_dict['INIT_INSTRUMENT_TYPE']
instruments_dict = instruments_dict.set_index('INDEX')
instruments_dict = instruments_dict.T.to_dict('dict')

global sub_adduct_pattern
sub_adduct_pattern = re.compile(r"\(|\)|(.*\[)|(\]([\d\+\-\*]*)?)")

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global ionization_mode_pattern
ionization_mode_pattern = re.compile(r"((?:^|\b)?APCI(?:\b|$)?)|((?:^|\b)?ACPI(?:\b|$)?)|((?:^|\b)?APPI(?:\b|$)?)|((?:^|\b)?EI(?:\b|$)?)|((?:^|\b)?ESI(?:\b|$)?)|((?:^|\b)?FAB(?:\b|$)?)|((?:^|\b)?MALDI(?:\b|$)?)",flags=re.IGNORECASE)

def normalize_empties(metadata_dict):
    """
    Normalizes empties in a metadata dictionary.

    :param metadata_dict: the dictionary containing metadata
    :return: the updated metadata dictionary
    """
    for k, v in metadata_dict.items():
        # For each item to replace
        if isinstance(v, str):
            if re.fullmatch(empty_pattern, v):
                metadata_dict[k] = ''

    return metadata_dict

def repair_inchi(metadata_dict):
    """
    :param metadata_dict: dictionary containing metadata
    :return: modified metadata dictionary with 'INCHI' key updated

    """
    inchi = metadata_dict['INCHI']

    if inchi:
        inchi = re.sub(repair_inchi_pattern,"InChI=",inchi)
        metadata_dict['INCHI'] = inchi

        return metadata_dict

    return metadata_dict

def repair_mol_descriptors(metadata_dict):
    """
        repair_mol_descriptors(metadata_dict)

        Repairs molecular descriptors in a dictionary containing SMILES, InChI, and InChIKey.
        Check if molecular descriptors are in the dedicated field.

        :param metadata_dict: A dictionary containing molecular descriptors.
        :return: The repaired dictionary with updated molecular descriptors.

        Example Usage:

        metadata_dict = {
            'SMILES': 'CC(=O)O',
            'INCHI': 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'
            'INCHIKEY': 'QTBSBXVTEAMEQO-UHFFFAOYSA-N'
        }

        repaired_dict = repair_mol_descriptors(metadata_dict)

        The repaired_dict will contain the repaired molecular descriptors.
    """
    smiles = metadata_dict['SMILES']
    inchi = metadata_dict['INCHI']
    inchikey = metadata_dict['INCHIKEY']

    if re.search(smiles_pattern, smiles) and re.search(inchi_pattern, inchi) and re.search(inchikey_pattern, inchikey):
        metadata_dict = repair_inchi(metadata_dict)
        return metadata_dict

    # SMILES
    if re.search(smiles_pattern, inchi):
        if not re.search(inchi_pattern, inchi) and not re.search(inchikey_pattern, inchi):
            metadata_dict['SMILES'] = inchi
            metadata_dict['INCHI'] = ''
    if re.search(smiles_pattern, inchikey):
        if not re.search(inchi_pattern, inchikey) and not re.search(inchikey_pattern, inchikey):
            metadata_dict['SMILES'] = inchikey
            metadata_dict['INCHIKEY'] = ''
    # INCHI
    if re.search(inchi_pattern, smiles):
        metadata_dict['INCHI'] = smiles
        metadata_dict['INCHI'] = ''
    if re.search(inchi_pattern, inchikey):
        metadata_dict['INCHI'] = inchikey
        metadata_dict['INCHIKEY'] = ''
    # INCHIKEY
    if re.search(inchikey_pattern, smiles):
        metadata_dict['INCHIKEY'] = smiles
        metadata_dict['SMILES'] = ''
    if re.search(inchikey_pattern, inchi):
        metadata_dict['INCHIKEY'] = inchi
        metadata_dict['INCHI'] = ''

    metadata_dict = repair_inchi(metadata_dict)

    return metadata_dict

def delete_no_smiles_no_inchi(metadata_dict):
    """
    Delete entries from the given metadata dictionary if both 'SMILES' and 'INCHI' keys have NaN values.

    :param metadata_dict: A dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary.
    :rtype: dict
    """
    if not metadata_dict["SMILES"] and not metadata_dict["INCHI"]:
        return None
    else:
        return metadata_dict

def normalize_adduct(metadata_dict):
    """
    Normalize adduct value in the given metadata dictionary.

    :param metadata_dict: The dictionary containing metadata information.
    :return: The modified metadata dictionary with normalized adduct value.
    """
    adduct = metadata_dict['PRECURSORTYPE']
    adduct = re.sub(sub_adduct_pattern, "", adduct)

    if adduct in adduct_dict:
        metadata_dict['PRECURSORTYPE'] = adduct_dict[adduct]

    return metadata_dict

def normalize_ionmode(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
        It should have the key "IONMODE" to represent the ionization mode.
    :return: The updated metadata dictionary with the "IONMODE" value normalized to either "positive" or "negative".
    """
    ionmode = metadata_dict["IONMODE"]

    if re.search(ionmode_pos_pattern, ionmode):
        ionmode = "positive"
    elif re.search(ionmode_neg_pattern, ionmode):
        ionmode = "negative"

    metadata_dict["IONMODE"] = ionmode

    return metadata_dict

def normalize_ms_level(metadata_dict):
    """
    Normalize the MS level value in the given metadata dictionary.

    :param metadata_dict: The dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized MS level value.
    :rtype: dict
    """
    ms_level = metadata_dict["MSLEVEL"]
    if ms_level:
        ms_level = re.findall(ms_level_pattern, ms_level)
        if ms_level:
            if len(ms_level) == 1:
                metadata_dict["MSLEVEL"] = ms_level[0]
            elif len(ms_level) >= 2:
                metadata_dict["MSLEVEL"] = f"{ms_level[0]}-{ms_level[1]}"
        else:
            metadata_dict["MSLEVEL"] = "2" # Init MSLEVEL to 2 by default
    else:
        metadata_dict["MSLEVEL"] = "2"  # Init MSLEVEL to 2 by default

    return metadata_dict

def in_filename(filename):
    """
    :param filename: The name of the file to be checked for the presence of the string "MSMS_Public" and matching a specific pattern.
    :return: Returns True if the filename does not contain "MSMS_Public" and matches a specific pattern defined by the In_Silico_pattern. Returns False otherwise.

    """
    if "MSMS_Public" not in filename:
        if re.search(In_Silico_pattern, filename):
            return True

    return False

def normalize_predicted(metadata_dict):
    """
    Normalize the predicted field in the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata information.
    :return: The updated metadata dictionary with the normalized predicted field.
    """
    comment_field = metadata_dict["COMMENT"]
    predicted = metadata_dict["PREDICTED"]
    filename = metadata_dict["FILENAME"]

    if re.search(In_Silico_pattern, comment_field) or predicted == "true" or in_filename(filename):
        metadata_dict["PREDICTED"] = "true"
        return metadata_dict
    else:
        metadata_dict["PREDICTED"] = "false"
        return metadata_dict

def normalize_retentiontime(metadata_dict):
    """
    Normalize the retention time value in the given metadata dictionary.

    :param metadata_dict: The dictionary containing the metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary with the normalized retention time.
    :rtype: dict
    """
    retientiontime = metadata_dict["RETENTIONTIME"]

    match = re.search(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(m|min|minute|minutes|s|sec|second|seconds|ms|millisecond|milliseconds)(?:\W)?", retientiontime, flags=re.IGNORECASE)

    if match:
        time = match.group(1)
        unit = match.group(2).lower()

        if not unit:
            retientiontime = str(float(time))
            metadata_dict["RETENTIONTIME"] = retientiontime
            return metadata_dict
        else:
            if unit in ["m", "min", "minute", "minutes"]:
                retientiontime = str(float(time))
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict
            elif unit in ["s", "sec", "second", "seconds"]:
                retientiontime = str(float(time) / 60)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict
            elif unit in ["ms", "millisecond", "milliseconds"]:
                retientiontime = str(float(time) / 60000)
                metadata_dict["RETENTIONTIME"] = retientiontime
                return metadata_dict

    return metadata_dict

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

def normalize_ionization(metadata_dict):
    """
    Normalize the ionization mode in the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata.
    :type metadata_dict: dict
    :return: The modified metadata dictionary.
    :rtype: dict
    """
    ionization_mode = re.search(ionization_mode_pattern, metadata_dict["IONIZATION"])

    if ionization_mode:
        ionization_mode = ionization_mode.group(1)
        if ionization_mode == "ACPI":
            ionization_mode = "APCI"
        metadata_dict["IONIZATION"] = ionization_mode
    else:
        ionization_mode_in_INSTRUMENTTYPE = re.search(ionization_mode_pattern, metadata_dict["INSTRUMENTTYPE"])
        if ionization_mode_in_INSTRUMENTTYPE:
            ionization_mode_in_INSTRUMENTTYPE = ionization_mode_in_INSTRUMENTTYPE.group(1)
            if ionization_mode_in_INSTRUMENTTYPE == "ACPI":
                ionization_mode_in_INSTRUMENTTYPE = "APCI"
            metadata_dict["IONIZATION"] = ionization_mode_in_INSTRUMENTTYPE

    return metadata_dict

def get_closest_match(instrument_name, instrument_list):
    """
    Find the closest match for the given instrument name in the given instrument list.

    :param instrument_name: The name of the instrument to find the closest match for.
    :param instrument_list: The list of instruments to search for a match in.

    :return: A tuple containing the closest match, if the similarity score is greater than or equal to 80. Otherwise, returns None.
    """
    closest_match = process.extractOne(instrument_name, instrument_list)
    if closest_match[1] < 80:  # si le score de similarité est inférieur à 80
        return None
    return closest_match[0]

def normalize_instruments_and_resolution(metadata_dict):
    """
    Normalize the instrument metadata in the given dictionary.

    :param metadata_dict: A dictionary containing the instrument metadata.
    :return: The normalized instrument metadata dictionary.
    """
    if metadata_dict["INSTRUMENT"]:
        metadata_dict_instrument = metadata_dict["INSTRUMENT"].lower()+metadata_dict["INSTRUMENTTYPE"].lower()
        closest_instrument = get_closest_match(metadata_dict_instrument, instruments_list)
        if closest_instrument:
            metadata_dict["INSTRUMENT"] = f'{instruments_dict[closest_instrument]["REF_INSTRUMENT"]}-{instruments_dict[closest_instrument]["REF_MODELE"]}'
            metadata_dict["INSTRUMENTTYPE"] = f'{instruments_dict[closest_instrument]["REF_SPECTRUM_TYPE"]}-{instruments_dict[closest_instrument]["REF_IONISATION"]}-{instruments_dict[closest_instrument]["REF_INSTRUMENT_TYPE"]}'
            metadata_dict["RESOLUTION"] = f'{instruments_dict[closest_instrument]["REF_RESOLUTION"]}'
            return metadata_dict

    return metadata_dict

def normalize_values(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    metadata_dict = normalize_empties(metadata_dict)

    metadata_dict = repair_mol_descriptors(metadata_dict)

    metadata_dict = delete_no_smiles_no_inchi(metadata_dict)

    if metadata_dict:
        metadata_dict = normalize_adduct(metadata_dict)
        metadata_dict = missing_precursormz_re_calculation(metadata_dict)
        metadata_dict = normalize_ionmode(metadata_dict)
        metadata_dict = normalize_ms_level(metadata_dict)
        metadata_dict = normalize_predicted(metadata_dict)
        metadata_dict = normalize_retentiontime(metadata_dict)
        metadata_dict = normalize_ionization(metadata_dict)
        metadata_dict = normalize_instruments_and_resolution(metadata_dict)

    return metadata_dict
