from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import pandas as pd
import warnings
import os

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

def complete_from_pubchem_datas(CONCATENATE_DF):
    # Créer une copie de concatenate_df pour éviter de modifier l'original
    concatenate_df_copy = CONCATENATE_DF.copy()

    # Joindre pubchem_df à concatenate_df sur la colonne 'INCHIKEY'
    enriched_df = concatenate_df_copy.merge(
        pubchem_datas,
        on='INCHIKEY',
        suffixes=('', '_pubchem'),
        how='left'
    )

    # Pour chaque colonne commune autre que 'INCHIKEY', remplacer les valeurs de concatenate_df par celles de pubchem_df
    columns_to_update = ['INCHI', 'SMILES', 'FORMULA', 'NAME', 'EXACTMASS', 'AVERAGEMASS']
    for col in columns_to_update:
        enriched_df[col] = enriched_df[col + '_pubchem'].combine_first(enriched_df[col])

    # Retirer les colonnes additionnelles de pubchem_df
    enriched_df.drop(columns=[col + '_pubchem' for col in columns_to_update], inplace=True)

    return enriched_df
