import numpy as np
import re

def normalize_empties(metadata_dict):
    regex_to_replace = [re.compile(re.escape(str(item)), re.I) for item in [0,
                                                                            0.0,
                                                                            "",
                                                                            "na",
                                                                            "n/a",
                                                                            "nan",
                                                                            "unknown",
                                                                            "unknown",
                                                                            "none",
                                                                            "?",
                                                                            "unk",
                                                                            np.nan]]

    for k, v in metadata_dict.items():
        # For each item to replace
        for regex in regex_to_replace:
            # If the current value matches, replace it
            if regex.fullmatch(str(v)):
                metadata_dict[k] = np.nan

    return metadata_dict

def determining_charge(adduct):
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

    matches = re.findall(r"((\+|\-|^)\d+)|(\+|\-)",adduct)

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

def normalize_adduct(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The updated metadata dictionary.

    This method takes a metadata dictionary as input and normalizes the adduct value in the dictionary. It extracts the adduct value from the "PRECURSORTYPE" key of the dictionary and checks
    * if it is in the correct format. If it is already in the correct format, the method does nothing and returns the original dictionary. If it is not in the correct format, the method
    * attempts to parse and normalize the adduct value.

    The adduct value should be in the format [molecule]charge, where molecule can be alphanumeric with optional +, -, (, and ) characters, and charge can be a numeric value followed by +
    *, -, or * characters. If the adduct value is not in this format, the method tries to extract the molecule part and the charge part separately.

    If the molecule part does not end with a charge, the method uses the `determining_charge` function to determine the charge value and appends it to the molecule part. If the molecule
    * part already contains a * character, it is preserved. The updated adduct value is then stored back in the dictionary under the "PRECURSORTYPE" key.

    If the molecule part already ends with a charge, the method checks if it also ends with a * character. If it does, the charge part is appended to the molecule part without any modifications
    *. If it does not end with a * character, the charge part is appended to the molecule part followed by a * character. The updated adduct value is stored back in the dictionary under
    * the "PRECURSORTYPE" key.

    Finally, the normalized metadata dictionary is returned.
    """
    adduct = metadata_dict["PRECURSORTYPE"]


    try:
        match = re.search(r"(\[([A-Za-z0-9\+\-\(\)]*)\]((?:[0-9]*)?[\+\-\*])*)(?:\/|$)?", adduct)
    except:
        return metadata_dict

    if match: # Si deja le format correct, on ne fait rien
        if not match.group(3): # si pas de charge a la fin
            charge = determining_charge(adduct)
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
        match = re.search("([A-Za-z0-9\+\-\(\)\*]*)", adduct)
        if match:
            if not re.search("(\d)?([\+\-\*])$", adduct): # si pas de charge a la fin
                charge = determining_charge(adduct)
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
                charge = re.search("(\d)?([\+\-\*])$", adduct)
                if charge:
                    metadata_dict["PRECURSORTYPE"] = "[" + match.group() + "]" + charge.group()
                    return metadata_dict
                else:
                    return metadata_dict

def normalize_values(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    metadata_dict = normalize_empties(metadata_dict)

    # metadata_dict = delete_no_smiles_inchi_inchikey(metadata_dict)
    
    metadata_dict = normalize_adduct(metadata_dict)
    # metadata_dict = normalize_ionmode(metadata_dict)
    # metadata_dict = normalize_retention_time(metadata_dict)
    # metadata_dict = normalize_ms_level(metadata_dict)
    # metadata_dict = normalize_synonymes(metadata_dict)
    # metadata_dict = normalize_formula(metadata_dict)
    # metadata_dict = normalize_predicted(metadata_dict)
    # metadata_dict = normalize_db_informations(metadata_dict)

    return metadata_dict