import globals_vars
import re

"""
si ionmode pos: mais adduit neg, on supprime (vis versa)
"""
def check_for_bad_adduct(metadata_dict):
    if metadata_dict['IONMODE'] == 'positive':
        if metadata_dict['PRECURSORTYPE'].endswith('-'):
            return None
        else:
            return metadata_dict
    elif metadata_dict['IONMODE'] == 'negative':
        if metadata_dict['PRECURSORTYPE'].startswith('+'):
            return None
        else:
            return metadata_dict

    return metadata_dict