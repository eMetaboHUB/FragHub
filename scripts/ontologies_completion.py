from tqdm import tqdm
import pandas as pd
import os

global ontologies_df
ontologies_df = pd.read_csv(os.path.abspath("../datas/ontologies_dict.csv"),sep=";", encoding="UTF-8") # Remplacez 'your_file.csv' par le chemin de votre fichier
ontologies_df = ontologies_df.fillna("UNKNOWN")

def ontologies_completion(spectrum_list):
    """
    ontologies_completion(spectrum_list)

    Add default values for 'CLASSYFIRE_CLASS' and 'NPCLASSIF_PATHWAY' columns, count unique INCHIKEYs, and update ontology information.

    Arguments:
    - spectrum_list: A DataFrame containing spectra information.

    Returns:
    - Updated spectrum_list with ontology information.

    Steps:
    1. Adds columns 'CLASSYFIRE_CLASS' and 'NPCLASSIF_PATHWAY' with default value 'UNKNOWN'.
    2. Merges spectrum_list with ontologies_df based on INCHIKEY.
    3. Updates the progress bar based on the merge progress.
    4. Returns the updated spectrum_list.
    """
    # Ajouter les colonnes 'CLASSYFIRE_CLASS' et 'NPCLASSIF_PATHWAY' avec des valeurs par défaut 'UNKNOWN'
    spectrum_list['CLASSYFIRE_CLASS'] = "UNKNOWN"
    spectrum_list['NPCLASSIF_PATHWAY'] = "UNKNOWN"

    # Compter le nombre de INCHIKEY uniques dans spectrum_list
    num_keys = spectrum_list['INCHIKEY'].nunique()

    # Initialise la barre de progression tqdm
    pbar = tqdm(total=num_keys, colour="green", unit=" key", desc="{:>70}".format("updating ontologies"))

    # Fusionner spectrum_list avec ontologies_df sur la colonne 'INCHIKEY'
    completed_df = pd.merge(
        spectrum_list,
        ontologies_df[['INCHIKEY', 'CLASSYFIRE_CLASS', 'NPCLASSIF_PATHWAY']],
        on='INCHIKEY',
        how='left'
    )

    # Vérifier et mettre à jour les colonnes avec les valeurs fusionnées `CLASSYFIRE_CLASS` et `NPCLASSIF_PATHWAY`
    completed_df['CLASSYFIRE_CLASS'] = completed_df['CLASSYFIRE_CLASS_y'].combine_first(
        completed_df['CLASSYFIRE_CLASS_x'])
    completed_df['NPCLASSIF_PATHWAY'] = completed_df['NPCLASSIF_PATHWAY_y'].combine_first(
        completed_df['NPCLASSIF_PATHWAY_x'])

    # Supprimer les colonnes temporaires
    completed_df.drop(
        columns=['CLASSYFIRE_CLASS_x', 'CLASSYFIRE_CLASS_y', 'NPCLASSIF_PATHWAY_x', 'NPCLASSIF_PATHWAY_y'],
        inplace=True)

    # Mettre à jour la barre de progression
    pbar.update(num_keys)
    pbar.close()

    return completed_df