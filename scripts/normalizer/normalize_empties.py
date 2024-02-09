import re

global empty_pattern
empty_pattern = re.compile(r"(^CCS:( .*)?)|(^\$:00in-source( .*)?)|(^0( .*)?)|(^0\.0( .*)?)|(^$)|(^na( .*)?)|(^n/a( .*)?)|(^nan( .*)?)|(^unknown( .*)?)|(^unknow( .*)?)|(^none( .*)?)|(^\?( .*)?)|(^unk( .*)?)|(^x( .*)?)", flags=re.IGNORECASE)

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