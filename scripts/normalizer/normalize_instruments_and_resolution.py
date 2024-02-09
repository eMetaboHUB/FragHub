from fuzzywuzzy import process
import pandas as pd
import itertools
import os

global instruments_list
instruments_list = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instrument_types_list = instruments_list['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_list = instruments_list['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_list = [' '.join(pair) for pair in itertools.zip_longest(instruments_list, instrument_types_list, fillvalue='')]

global instruments_dict
instruments_dict = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instruments_dict['INIT_INSTRUMENT'] = instruments_dict['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_dict['INIT_INSTRUMENT_TYPE'] = instruments_dict['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_dict['INDEX'] = instruments_dict['INIT_INSTRUMENT'] + " " + instruments_dict['INIT_INSTRUMENT_TYPE']
instruments_dict = instruments_dict.set_index('INDEX')
instruments_dict = instruments_dict.T.to_dict('dict')

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
        metadata_dict_instrument = metadata_dict["INSTRUMENT"].lower() + " " + metadata_dict["INSTRUMENTTYPE"].lower()
        closest_instrument = get_closest_match(metadata_dict_instrument, instruments_list)
        if closest_instrument:
            metadata_dict["INSTRUMENT"] = f'{instruments_dict[closest_instrument]["REF_INSTRUMENT"]}-{instruments_dict[closest_instrument]["REF_MODELE"]}'
            metadata_dict["INSTRUMENTTYPE"] = f'{instruments_dict[closest_instrument]["REF_SPECTRUM_TYPE"]}-{instruments_dict[closest_instrument]["REF_IONISATION"]}-{instruments_dict[closest_instrument]["REF_INSTRUMENT_TYPE"]}'
            metadata_dict["RESOLUTION"] = f'{instruments_dict[closest_instrument]["REF_RESOLUTION"]}'
            return metadata_dict

    return metadata_dict