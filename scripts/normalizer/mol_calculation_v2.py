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