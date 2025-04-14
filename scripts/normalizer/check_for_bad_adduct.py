import scripts.deletion_report
import scripts.globals_vars
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
    ionmode = metadata_dict['IONMODE']

    instrument_type = metadata_dict["INSTRUMENTTYPE"]
    if re.search(r"\bEI\b", instrument_type):
        if not adduct:
            return metadata_dict

    if adduct == "M":
        if ionmode == 'positive':
            adduct = "[M]+"
            metadata_dict['PRECURSORTYPE'] = adduct
            return metadata_dict
        elif ionmode == 'negative':
            adduct = "[M]-"
            metadata_dict['PRECURSORTYPE'] = adduct
            return metadata_dict


    if not re.search(globals_vars.is_adduct_pattern, adduct):
        metadata_dict['DELETION_REASON'] = "spectrum deleted because its adduct field is empty or the value entered is not an adduct"
        deletion_report.deleted_spectrum_list.append(metadata_dict)
        deletion_report.no_or_bad_adduct += 1
        return None

    if ionmode == 'positive':
        if adduct in globals_vars.adduct_massdiff_dict_NEG:
            metadata_dict['DELETION_REASON'] = "spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode)."
            deletion_report.deleted_spectrum_list.append(metadata_dict)
            deletion_report.no_or_bad_adduct += 1
            return None
        else:
            return metadata_dict
    elif ionmode == 'negative':
        if adduct in globals_vars.adduct_massdiff_dict_POS:
            metadata_dict['DELETION_REASON'] = "spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode)."
            deletion_report.deleted_spectrum_list.append(metadata_dict)
            deletion_report.no_or_bad_adduct += 1
            return None
        else:
            return metadata_dict

    return metadata_dict