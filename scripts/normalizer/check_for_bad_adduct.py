import deletion_report
import globals_vars
import re

def check_for_bad_adduct(metadata_dict):
    """
    Check for inconsistent adducts based on ionization mode.

    This function validates the consistency of the ionization mode and precursor type
    defined in a given metadata dictionary. It ensures that positive ion modes do not
    contain precursor types ending with a negative sign and that negative ion modes do
    not contain precursor types starting with a positive sign. If the adduct is found
    to be inconsistent, the function returns None. If there are no issues, it returns
    the original metadata dictionary.

    Arguments:
        metadata_dict (dict): A dictionary containing ionization mode and precursor
        type information. The keys 'IONMODE' and 'PRECURSORTYPE' should be defined.

    Returns:
        dict or None: The same metadata dictionary if no inconsistencies are found;
        otherwise, None.
    """
    adduct = metadata_dict['PRECURSORTYPE']

    if not re.search(globals_vars.is_adduct_pattern, adduct):
        deletion_report.no_or_bad_adduct += 1
        return None

    if metadata_dict['IONMODE'] == 'positive':
        if metadata_dict['PRECURSORTYPE'].endswith('-'):
            deletion_report.no_or_bad_adduct += 1
            return None
        else:
            return metadata_dict
    elif metadata_dict['IONMODE'] == 'negative':
        if metadata_dict['PRECURSORTYPE'].endswith('+'):
            deletion_report.no_or_bad_adduct += 1
            return None
        else:
            return metadata_dict

    return metadata_dict