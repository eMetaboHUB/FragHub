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

def correct_and_convert_inchi_smiles(row):
    # Si INCHIKEY est présent et correspond au pattern, ne pas faire d'autres traitements
    if pd.notna(row['INCHIKEY']) and inchikey_pattern.match(row['INCHIKEY']):
        return row

    if pd.notna(row['INCHI']):
        if 'InChI=' not in row['INCHI']:
            cleaned_value = re.sub(indigo_smiles_correction_pattern, "", row['INCHI'])
            if Chem.MolFromSmiles(cleaned_value):
                row['INCHI'] = Chem.MolToSmiles(Chem.MolFromSmiles(cleaned_value))
            else:
                row['INCHI'] = ''  # Nettoyer si ce n'est ni un INCHI valide, ni un SMILES valide

    if pd.notna(row['SMILES']):
        if 'InChI=' in row['SMILES']:
            if Chem.MolFromInchi(row['SMILES']):
                row['SMILES'] = Chem.MolToInchi(Chem.MolFromInchi(row['SMILES']))
            else:
                row['SMILES'] = ''  # Nettoyer si ce n'est pas un SMILES valide

    return row

def mols_derivation_and_calculation(CONCATENATE_DF):
    # Application de la fonction à chaque ligne du DataFrame
    CONCATENATE_DF = CONCATENATE_DF.apply(correct_and_convert_inchi_smiles, axis=1)