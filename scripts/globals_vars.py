import pandas as pd
import psutil
import os
import re


global available_memory
available_memory = psutil.virtual_memory().available

global cpu_count
cpu_count = os.cpu_count()  # Nombre de cœurs logiques

global peak_list_csv_to_json_pattern
peak_list_csv_to_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global peak_list_json_to_json_pattern
peak_list_json_to_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global keys_dict
Key_dataframe = pd.read_csv(os.path.abspath("../datas/key_to_convert.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
keys_dict = dict(zip(Key_dataframe['known_synonym'], Key_dataframe['fraghub_default'].str.upper()))

global keys_list
keys_list = ['FILENAME',
             'PREDICTED',
             'SPLASH',
             'SPECTRUMID',
             'RESOLUTION',
             'SYNON',
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
             'ENTROPY',
             'CLASSYFIRE_SUPERCLASS',
             'CLASSYFIRE_CLASS',
             'CLASSYFIRE_SUBCLASS',
             'NPCLASS_PATHWAY',
             'NPCLASS_SUPERCLASS',
             'NPCLASS_CLASS',
             'NUM PEAKS',
             'PEAKS_LIST']