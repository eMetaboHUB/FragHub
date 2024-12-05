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
    CONCATENATE_DF['calculation'] = ''
    # Ajouter une colonne 'INDEX' qui numérote les lignes
    CONCATENATE_DF['INDEX'] = range(len(CONCATENATE_DF))

    # Supprimer temporairement les avertissements
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)

        def update_from_inchikey(df1, df2):
            # Barre de progression "fake"
            progress_bar = tqdm(df2, unit=" rows", colour="green", desc="{:>70}".format("completting from pubchem datas"))
            progress_bar.update(1)  # Passe directement à 100%
            progress_bar.close()

            # Obtenez l'intersection des INCHIKEY présents dans les deux DataFrames
            common_inchikeys = set(df1['INCHIKEY']).intersection(df2['INCHIKEY'])
            # Créer un nouvel ensemble sans éléments vides ou NaN
            common_inchikeys = {inchikey for inchikey in common_inchikeys if inchikey and pd.notna(inchikey)}

            # Filtrer les DataFrames pour ne conserver que les INCHIKEYs communs
            df1_filtered = df1[df1['INCHIKEY'].isin(common_inchikeys)].copy()
            df2_filtered = df2[df2['INCHIKEY'].isin(common_inchikeys)].copy()

            df1_filtered.set_index('INCHIKEY', inplace=True)
            df2_filtered.set_index('INCHIKEY', inplace=True)

            # Mettre à jour df1_filtered dans df2_filtered
            df2_filtered.update(df1_filtered, errors='ignore')
            df2_filtered.reset_index(inplace=True)

            # Assurez-vous de restaurer l'index original pour df2 avant d'effectuer la mise à jour finale
            df2.set_index('INDEX', inplace=True)

            # Mettre à jour uniquement les parties de df2_filtered qui ont l'index correspondant
            df2.update(df2_filtered.set_index('INDEX'))

            # Réinitialisez l'index pour que 'INDEX' soit une colonne à nouveau
            df2.reset_index(inplace=True)

        # Effectuer d'abord la mise à jour par 'INCHIKEY'
        update_from_inchikey(pubchem_datas, CONCATENATE_DF)

    return CONCATENATE_DF
