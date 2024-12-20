from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import psutil
import json
import os
import re

# =================================================== REGEX PATTERN ====================================================

# ============ Parsors regex ============
global metadata_strip_value_pattern
metadata_strip_value_pattern = re.compile(r"^\"|\"$")

global metadata_fields_name_pattern
metadata_fields_name_pattern = re.compile(r'^[\W_]+|[\W_]+$')

global metadata_pattern_mgf
metadata_pattern_mgf = re.compile(r"([^:\n]*?)=\s*([^\n]*)(?:\n|$)")

global metadata_pattern_msp
metadata_pattern_msp = re.compile(r"([^:]*):(?: )?([^\n]*)(?:\n|$)")

global computed_pattern
computed_pattern = re.compile(r"computed", flags=re.IGNORECASE)

global comment_pattern
comment_pattern = re.compile(r'comment.*', flags=re.IGNORECASE)

global peak_list_split_pattern
peak_list_split_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global peak_list_json_pattern
peak_list_json_pattern = re.compile(r"(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)(?:\s+|:|,|, )(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global sub_fields_pattern
sub_fields_pattern = re.compile(r"(\S+?)=\"([^\"]*)\"|\"(\w+?)=([^\"]*)\"|\"([^\"]*?)=([^\"]*)\"|(\S+?)=(\d+(?:[.,]\d*)?)|(\S+?)=(.*?)(?:;|\n|$)")

global metadata_peak_list_split_pattern_mgf
metadata_peak_list_split_pattern_mgf = re.compile(r"([\s\S]*=.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")

global metadata_peak_list_split_pattern_msp
metadata_peak_list_split_pattern_msp = re.compile(r"([\s\S]*:.*[0-9]*\n)(((-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(\s+|:)(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(.*)(\n|$))*)")
# ======================================

# ===== normalizers regex pattern ======
global indigo_smiles_correction_pattern
indigo_smiles_correction_pattern = re.compile(r"\|[\s\S]*")

global sub_adduct_pattern
sub_adduct_pattern = re.compile(r"\(|\)|(.*\[)|(\]([\d\+\-\*]*)?)")

global float_check_pattern
float_check_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)")

global ms_level_pattern
ms_level_pattern = re.compile(r"(?:ms)?(\d)", flags=re.IGNORECASE)

global ionmode_pos_pattern
ionmode_pos_pattern = re.compile(r"^p|^\+|^pos", flags=re.IGNORECASE)

global ionmode_neg_pattern
ionmode_neg_pattern = re.compile(r"^n|^\-|^neg", flags=re.IGNORECASE)

global repair_inchi_pattern
repair_inchi_pattern = re.compile(r"^(inchi=)?", flags=re.IGNORECASE)

global inchi_pattern
inchi_pattern = re.compile(r"InChI=.*|\/[0-9A-Z]*\/", flags=re.IGNORECASE) # Match inchi

global smiles_pattern
smiles_pattern = re.compile(r"[^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,}", flags=re.IGNORECASE) # Match smiles

global inchikey_pattern
inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey

global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)

global retention_time_pattern
retention_time_pattern = re.compile(r"(-?\d+[.,]?\d*(?:[Ee][+-]?\d+)?)(?:\W)?(m|min|minute|minutes|s|sec|second|seconds|ms|millisecond|milliseconds)(?:\W)?")

global ionization_mode_pattern
ionization_mode_pattern = re.compile(r"((?:^|\b)?APCI(?:\b|$)?)|((?:^|\b)?ACPI(?:\b|$)?)|((?:^|\b)?APPI(?:\b|$)?)|((?:^|\b)?EI(?:\b|$)?)|((?:^|\b)?ESI(?:\b|$)?)|((?:^|\b)?FAB(?:\b|$)?)|((?:^|\b)?MALDI(?:\b|$)?)",flags=re.IGNORECASE)

global empty_pattern
empty_pattern = re.compile(r"(^CCS:( .*)?)|(^\$:00in-source( .*)?)|(^0( .*)?)|(^0\.0( .*)?)|(^$)|(^na( .*)?)|(^n/a( .*)?)|(^nan( .*)?)|(^unknown( .*)?)|(^unknow( .*)?)|(^none( .*)?)|(^\?( .*)?)|(^unk( .*)?)|(^x( .*)?)", flags=re.IGNORECASE)
# =====================================

# ======================================================================================================================

# =================================================== READ FILES =======================================================

global ontologies_df
files = [f for f in os.listdir(os.path.abspath("../datas/ontologies_datas")) if 'ontologies_dict' in f]
ontologies_df = pd.concat((pd.read_csv(os.path.join(os.path.abspath("../datas/ontologies_datas/"), f), sep=";", encoding="UTF-8") for f in files), ignore_index=True)

# ================

# Dossier contenant les fichiers CSV
folder_path = '../datas/pubchem_datas/'
# Liste pour stocker chaque DataFrame
all_dfs = []
# Fonction pour lire un fichier CSV
def read_csv(file_path):
    return pd.read_csv(file_path, sep=';', quotechar='"', encoding='utf-8')

# Utiliser ThreadPoolExecutor pour lire les fichiers en parallèle
with ThreadPoolExecutor() as executor:
    futures = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            futures.append(executor.submit(read_csv, file_path))
    # Récupérer les résultats des futures
    for future in futures:
        all_dfs.append(future.result())

# Concaténer tous les DataFrames
global pubchem_datas
pubchem_datas = pd.concat(all_dfs, ignore_index=True)

# ================

global adduct_dict, adduct_massdiff_dict
adduct_dataframe = pd.read_csv(os.path.abspath("../datas/adduct_to_convert.csv"), sep=";", encoding="UTF-8")
adduct_dict = dict(zip(adduct_dataframe['known_adduct'], adduct_dataframe['fraghub_default']))
adduct_massdiff_dict = dict(zip(adduct_dataframe['fraghub_default'], adduct_dataframe['massdiff']))

# ================

global instrument_tree
with open('../datas/instruments_tree.json', 'r') as f:
    instrument_tree = json.load(f)

# ================

global keys_dict
Key_dataframe = pd.read_csv(os.path.abspath("../datas/key_to_convert.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
keys_dict = dict(zip(Key_dataframe['known_synonym'], Key_dataframe['fraghub_default'].str.upper()))

# ======================================================================================================================

# =====================================================LIST=============================================================
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

# ======================================================================================================================

# ====================================================OTHER VARS========================================================

global available_memory
available_memory = psutil.virtual_memory().available

global cpu_count
cpu_count = os.cpu_count()  # Nombre de cœurs logiques

# ======================================================================================================================
