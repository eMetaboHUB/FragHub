import numpy as np
import re

global repair_inchi_pattern
repair_inchi_pattern = re.compile("^(inchi=)?", flags=re.IGNORECASE)

global smiles_pattern
smiles_pattern = re.compile("[^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,}", flags=re.IGNORECASE) # Match smiles

global inchi_pattern
inchi_pattern = re.compile("InChI=.*|\/[0-9A-Z]*\/", flags=re.IGNORECASE) # Match inchi

global inchikey_pattern
inchikey_pattern = re.compile("([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey

global adduct_pattern
adduct_pattern = re.compile(r"(\[([A-Za-z0-9\+\-\(\)]*)\]((?:[0-9]*)?[\+\-\*])*)(?:\/|$)?")

global charge_pattern
charge_pattern = re.compile(r"((\+|\-|^)\d+)|(\+|\-)")

global adduct_pattern_2
adduct_pattern_2 = re.compile("([A-Za-z0-9\+\-\(\)\*]*)")

global ending_by_charge_pattern
ending_by_charge_pattern = re.compile("(\d)?([\+\-\*])$")

global ionmode_pos_pattern
ionmode_pos_pattern = re.compile("^p|^\+|^pos", flags=re.IGNORECASE)

global ionmode_neg_pattern
ionmode_neg_pattern = re.compile("^n|^\-|^neg", flags=re.IGNORECASE)

global precursortype_pos_pattern
precursortype_pos_pattern = re.compile("\][\+\*]*$")

global precursortype_neg_pattern
precursortype_neg_pattern = re.compile("\][\-\*]*$")

global ms_level_pattern
ms_level_pattern = re.compile("(?:ms)?(\d)", flags=re.IGNORECASE)


def normalize_empties(metadata_dict):
    """
    Normalizes empties in a metadata dictionary.

    :param metadata_dict: the dictionary containing metadata
    :return: the updated metadata dictionary
    """
    regex_to_replace = [re.compile(re.escape(str(item)), re.I) for item in [0, 0.0, "", "na", "n/a", "nan", "unknown", "unknow", "none", "?", "unk", np.nan]]

    for k, v in metadata_dict.items():
        # For each item to replace
        for regex in regex_to_replace:
            # If the current value matches, replace it
            if regex.fullmatch(str(v)):
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




def determining_adduct_charge(adduct): # NOTE: ATTENTION: pas certains que ce soit correcte
    """
    Calculate the charge of a given adduct.

    :param adduct: The adduct string to calculate the charge for.
                   The adduct string should contain a combination of numbers and
                   '+' or '-' symbols to represent the charge. For example, '+2',
                   '-1', '+', '-'.
    :return: Return the calculated charge as a string. The returned string represents
             the charge and can have one of the following formats: '+', '-', '+X',
             '-X', where X is a positive integer.
    """
    sum = 0
    min = 0
    max = 0

    matches = re.findall(charge_pattern,adduct)

    # determine signe
    if matches:
        for match in matches:
            charge = next((item for item in match if item), '')
            if charge == "-":
                sum -= 1
                if min > -1:
                    min = -1
            elif charge == "+":
                sum += 1
                if max < 1:
                    max = 1
            else:
                sum += int(charge)
                if int(charge) > max:
                    max = int(charge)
                if int(charge) < min:
                    min = int(charge)

        # determine charge
        if sum > 0:
            if max == 1:
                return "+"
            else:
                return str(max)+"+"
        elif sum < 0:
            if min == -1:
                return "-"
            else:
                return str(abs(min))+"-"
    else:
        return None


def delete_no_smiles_inchi_inchikey(metadata_dict):
    """
    Delete entries from the given metadata dictionary if both 'SMILES' and 'INCHI' keys have NaN values.

    :param metadata_dict: A dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary.
    :rtype: dict
    """
    if not metadata_dict.get('SMILES') and not metadata_dict.get('INCHI'):
        return None
    else:
        return metadata_dict



def normalize_adduct(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The modified metadata dictionary.

    This method takes a dictionary of metadata information and normalizes the "PRECURSORTYPE" value according to a specific format. The method checks if the value already matches the format
    *, and if not, it modifies it accordingly.

    If the value already matches the expected format, the method does nothing and returns the unmodified metadata dictionary. If the value does not match the format, the method attempts
    * to determine the charge value and add it to the end of the value in the correct format.

    Note: The method relies on an external function called determining_charge() to determine the charge value. This function is not included in this documentation.

    Example usage:
    metadata_dict = {"PRECURSORTYPE": "[M+H]"}
    normalized_dict = normalize_adduct(metadata_dict)
    print(normalized_dict)

    Output:
    {"PRECURSORTYPE": "[M+H]"}
    """
    adduct = metadata_dict["PRECURSORTYPE"]


    match = re.search(adduct_pattern, adduct)

    if match: # Si deja le format correct, on ne fait rien
        if not match.group(3): # si pas de charge a la fin
            charge = determining_adduct_charge(adduct)
            if "*" in match.group():
                if charge:
                    metadata_dict["PRECURSORTYPE"] = adduct + charge + "*"
                    return metadata_dict
                else:
                    return metadata_dict
            else:
                if charge:
                    metadata_dict["PRECURSORTYPE"] = adduct + charge
                    return metadata_dict
                else:
                    return metadata_dict
        return metadata_dict
    else: # pas le format correct
        match = re.search(adduct_pattern_2, adduct)
        if match:
            if not re.search(ending_by_charge_pattern, adduct): # si pas de charge a la fin
                charge = determining_adduct_charge(adduct)
                if "*" in match.group():
                    if charge:
                        metadata_dict["PRECURSORTYPE"] = "[" + match.group() + "]" + charge + "*"
                        return metadata_dict
                    else:
                        return metadata_dict
                else:
                    if charge:
                        metadata_dict["PRECURSORTYPE"] = "[" + match.group() + "]" + charge
                        return metadata_dict
                    else:
                        return metadata_dict
            else:
                charge = re.search(ending_by_charge_pattern, adduct)
                if charge:
                    metadata_dict["PRECURSORTYPE"] = "[" + match.group() + "]" + charge.group()
                    return metadata_dict
                else:
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
        ms_level = re.search(ms_level_pattern, ms_level)
        if ms_level:
            metadata_dict["MSLEVEL"] = ms_level.group(1)

    return metadata_dict

def normalize_values(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    metadata_dict = normalize_empties(metadata_dict)

    metadata_dict = repair_mol_descriptors(metadata_dict)

    metadata_dict = delete_no_smiles_inchi_inchikey(metadata_dict)

    if metadata_dict:
         metadata_dict = normalize_adduct(metadata_dict)
         metadata_dict = normalize_ionmode(metadata_dict)
         # metadata_dict = normalize_retention_time(metadata_dict)
         metadata_dict = normalize_ms_level(metadata_dict)
         # metadata_dict = normalize_synonymes(metadata_dict)
         # metadata_dict = normalize_formula(metadata_dict)
         # metadata_dict = normalize_predicted(metadata_dict)
         # metadata_dict = normalize_db_informations(metadata_dict)

    return metadata_dict
