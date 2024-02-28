import json
import re

global instrument_tree
with open('../datas/instruments_tree.json', 'r') as f:
    instrument_tree = json.load(f)

def clean_instrument(instrument):
    """
    Clean the instrument name string by removing specific prefixes and modifying the format.

    :param instrument: The instrument name string to be cleaned.
    :return: The cleaned instrument name string.
    """
    instrument = re.sub("-tof", "tof", instrument)
    instrument = re.sub("q-", "q", instrument)
    instrument = re.sub("q exactive", "qexactive", instrument)
    instrument = re.sub("applied biosystems", "ab sciex", instrument)
    instrument = re.sub(" ab ", " ab sciex ", instrument)
    instrument = re.sub("sciex", " ab sciex ", instrument)
    instrument = re.sub("triple(-| )?tof", "qqq", instrument)
    instrument = re.sub("triple(-| )?quad", "qqq", instrument)
    instrument = re.sub("... uplc ...", "", instrument)

    return instrument

def clean_instrument_type(instrument_type):
    """
    Cleans the given instrument string by removing specific substrings and replacing hyphens with spaces.

    :param instrument_str: The instrument string to be cleaned.
    :type instrument_str: str
    :return: The cleaned instrument string.
    :rtype: str
    """
    instrument_type = re.sub("-tof", "tof", instrument_type)
    instrument_type = re.sub("q-", "q", instrument_type)
    instrument_type = re.sub("-", " ", instrument_type)
    instrument_type = re.sub("q exactive", "qexactive", instrument_type)
    instrument_type = re.sub("applied biosystems", "ab sciex", instrument_type)
    instrument_type = re.sub(" ab ", " ab sciex ", instrument_type)
    instrument_type = re.sub("sciex", " ab sciex ", instrument_type)
    instrument_type = re.sub("triple(-| )?tof", "qqq", instrument_type)
    instrument_type = re.sub("triple(-| )?quad", "qqq", instrument_type)
    instrument_type = re.sub("... uplc ...", "", instrument_type)

    return instrument_type

def clean_spectrum_instrument_info(metadata_dict):
    """
    Cleans the spectrum instrument information from the given metadata dictionary.

    :param metadata_dict: The dictionary containing the metadata information.
    :return: The cleaned instrument information.
    """
    instrument = metadata_dict['INSTRUMENT'].lower()
    instrument_type = metadata_dict["INSTRUMENTTYPE"].lower()

    instrument = clean_instrument(instrument)
    instrument_type = clean_instrument_type(instrument_type)


    instrument_infos = instrument + " " + instrument_type
    instrument_infos = re.sub(r'[^-\w\s]', ' ', instrument_infos)
    instrument_infos = ' '.join(instrument_infos.split()).strip()

    return instrument_infos

def normalize_instruments_and_resolution(metadata_dict):
    """
    Normalize the instrument metadata in the given dictionary.

    :param metadata_dict: A dictionary containing the instrument metadata.
    :return: The normalized instrument metadata dictionary.
    """
    clean_spectrum_instrument_info(metadata_dict)

    return metadata_dict