import pandas as pd
import os

global keys_dict
Key_dataframe = pd.read_csv(os.path.abspath("../datas/key_to_convert.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
keys_dict = dict(zip(Key_dataframe['known_synonym'], Key_dataframe['fraghub_default'].str.upper()))

global keys_list
keys_list = ['FILENAME',
             'PREDICTED',
             'FRAGHUBID',
             'SPECTRUMID',
             'RESOLUTION',
             'SYNON',
             'CHARGE',
             'IONIZATION',
             'MSLEVEL',
             'FRAGMENTATIONMODE',
             'NAME',
             'PRECURSORMZ',
             'EXACTMASS',
             'AVERAGEMASS',
             'PRECURSORTYPE',
             'INSTRUMENTTYPE',
             'INSTRUMENT',
             'SMILES',
             'INCHI',
             'INCHIKEY',
             'COLLISIONENERGY',
             'FORMULA',
             'RETENTIONTIME',
             'IONMODE',
             'COMMENT',
             'NUM PEAKS',
             'PEAKS_LIST']

def convert_keys(metadata_dict):
    """
    Convert keys in metadata_dict based on the provided keys_dict and keys_list.

    :param metadata_dict: A dictionary containing metadata information.
    :return: A dictionary with converted keys based on the provided keys_dict and keys_list.
    """
    converted = {keys_dict[key.lower()]: val for key, val in metadata_dict.items() if key.lower() in keys_dict and keys_dict[key.lower()] in keys_list}

    converted.update({key: "" for key in keys_list if key not in converted})

    return converted