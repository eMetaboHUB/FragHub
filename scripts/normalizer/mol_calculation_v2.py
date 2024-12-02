from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from concurrent.futures import ThreadPoolExecutor
from rdkit import RDLogger, Chem
from tqdm import tqdm
import pandas as pd
import re
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

# Dossier contenant les fichiers CSV
folder_path = '../datas/pubchem_datas'
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
pubchem_datas = pd.concat(all_dfs, ignore_index=True)

inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey
indigo_smiles_correction_pattern = re.compile(r"\|[\s\S]*")

def merge_datas_from_pubchem(CONCATENATE_DF, pubchem_datas):
    # Premier merge par INCHIKEY
    CONCATENATE_DF = pd.merge(CONCATENATE_DF, pubchem_datas, on='INCHIKEY', how='left', suffixes=('', '_pubchem'))

    # Supprimer les colonnes suffixées après l'assignation
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    # Masque pour les lignes où 'calculation' n'est pas 'done'
    not_done_mask = CONCATENATE_DF['calculation'] != 'done'

    # Deuxième fusion par INCHI sur les lignes non encore marquées comme 'done'
    CONCATENATE_DF.loc[not_done_mask] = pd.merge(
        CONCATENATE_DF[not_done_mask],
        pubchem_datas,
        on='INCHI',
        how='left',
        suffixes=('', '_pubchem')
    )

    # Supprimer à nouveau les colonnes suffixées après l'assignation
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    # Mettre à jour le masque pour les lignes restantes
    not_done_mask = CONCATENATE_DF['calculation'] != 'done'

    # Troisième fusion par SMILES sur les lignes non encore marquées comme 'done'
    CONCATENATE_DF.loc[not_done_mask] = pd.merge(
        CONCATENATE_DF[not_done_mask],
        pubchem_datas,
        on='SMILES',
        how='left',
        suffixes=('', '_pubchem')
    )

    # Supprimer les colonnes suffixées après la troisième fusion
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    return CONCATENATE_DF

def mols_derivation_and_calculation(CONCATENATE_DF):
    CONCATENATE_DF['calculation'] = ''

    CONCATENATE_DF = merge_datas_from_pubchem(CONCATENATE_DF, pubchem_datas)

