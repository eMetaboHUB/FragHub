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
    2. Counts the number of unique INCHIKEYs in spectrum_list.
    3. Initializes a tqdm progress bar for tracking updates.
    4. Defines and applies fill_row function to update ontology information based on INCHIKEY matches.
    5. Closes the progress bar and returns the updated spectrum_list.
    """
    # Ajouter les colonnes 'CLASSYFIRE_CLASS' et 'NPCLASSIF_PATHWAY' avec des valeurs par défaut 'UNKNOWN'
    spectrum_list['CLASSYFIRE_CLASS'] = "UNKNOWN"
    spectrum_list['NPCLASSIF_PATHWAY'] = "UNKNOWN"

    # Compter le nombre de INCHIKEY uniques dans spectrum_list
    num_keys = spectrum_list['INCHIKEY'].nunique()

    # Initialise la barre de progression tqdm
    pbar = tqdm(total=num_keys, colour="green", unit=" key", desc="{:>70}".format("updating ontologies"))

    def fill_row(row):
        # Mise à jour des valeurs si les INCHIKEY correspondent
        match = ontologies_df[ontologies_df['INCHIKEY'] == row['INCHIKEY']]
        if not match.empty:
            row['CLASSYFIRE_CLASS'] = match['CLASSYFIRE_CLASS'].values[0]
            row['NPCLASSIF_PATHWAY'] = match['NPCLASSIF_PATHWAY'].values[0]
        pbar.update()  # Mise à jour de la barre de progression
        return row

    # Appliquer la fonction à chaque ligne de spectrum_list
    spectrum_list = spectrum_list.apply(fill_row, axis=1)

    pbar.close()  # Fermer la barre de progression

    return spectrum_list