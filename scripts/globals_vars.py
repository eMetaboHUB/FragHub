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

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*=.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:\n]*?)=\s*([^\n]*)(?:\n|$)")

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global metadata_peak_list_split_pattern
metadata_peak_list_split_pattern = re.compile(r"([\s\S]*:.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_pattern
metadata_pattern = re.compile(r"([^:]*):(?: )?([^\n]*)(?:\n|$)")

global computed_pattern
computed_pattern = re.compile(r"computed", flags=re.IGNORECASE)

global comment_pattern
comment_pattern = re.compile(r'comment.*', flags=re.IGNORECASE)

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global sub_fields_pattern
sub_fields_pattern = re.compile(r"(\S+?)=\"([^\"]*)\"|\"(\w+?)=([^\"]*)\"|\"([^\"]*?)=([^\"]*)\"|(\S+?)=(\d+(?:[.,]\d*)?)|(\S+?)=(.*?)(?:;|\n|$)")

global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

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

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global adduct_massdiff_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_massdiff_dict = dict(zip(adduct_dataframe['fraghub_default'], adduct_dataframe['massdiff']))