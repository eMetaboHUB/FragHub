import pandas as pd
import itertools
import os
import re

global instruments_list
instruments_list = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instrument_types_list = instruments_list['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_list = instruments_list['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_list = [' '.join(re.sub(r'[^\w\s]', ' ', ' '.join(pair)).split()) for pair in itertools.zip_longest(instruments_list, instrument_types_list, fillvalue='')]

global instruments_dict
instruments_dict = pd.read_csv(os.path.abspath("../datas/instruments_catalogue.csv"), sep=";", encoding="UTF-8")
instruments_dict['INIT_INSTRUMENT'] = instruments_dict['INIT_INSTRUMENT'].str.lower().fillna('')
instruments_dict['INIT_INSTRUMENT_TYPE'] = instruments_dict['INIT_INSTRUMENT_TYPE'].str.lower().fillna('')
instruments_dict['INDEX'] = instruments_dict.apply(lambda row: ' '.join([''.join(sorted(word)) for word in re.sub(r'[^\w\s]', ' ', row['INIT_INSTRUMENT'] + " " + row['INIT_INSTRUMENT_TYPE']).split()]), axis=1)
instruments_dict = instruments_dict.set_index('INDEX')
instruments_dict = instruments_dict.T.to_dict('dict')

def normalize_instruments_and_resolution(metadata_dict):
    """
    Normalize the instrument metadata in the given dictionary.

    :param metadata_dict: A dictionary containing the instrument metadata.
    :return: The normalized instrument metadata dictionary.
    """
    if metadata_dict["INSTRUMENT"]:
        metadata_dict_instrument = ' '.join(re.sub(r'[^\w\s]', ' ', metadata_dict["INSTRUMENT"].lower() + " " + metadata_dict["INSTRUMENTTYPE"].lower()).split())
        normalized_instrument_name = ''.join(sorted(list(metadata_dict_instrument)))

        if normalized_instrument_name in instruments_dict:
            metadata_dict["INSTRUMENT"] = f'{instruments_dict[normalized_instrument_name]["REF_INSTRUMENT"]}-{instruments_dict[normalized_instrument_name]["REF_MODELE"]}'
            metadata_dict["INSTRUMENTTYPE"] = f'{instruments_dict[normalized_instrument_name]["REF_SPECTRUM_TYPE"]}-{instruments_dict[normalized_instrument_name]["REF_IONISATION"]}-{instruments_dict[normalized_instrument_name]["REF_INSTRUMENT_TYPE"]}'
            metadata_dict["RESOLUTION"] = f'{instruments_dict[normalized_instrument_name]["REF_RESOLUTION"]}'
            return metadata_dict

    return metadata_dict