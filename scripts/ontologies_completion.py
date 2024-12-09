from tqdm import tqdm
import pandas as pd
import os

global ontologies_df
files = [f for f in os.listdir(os.path.abspath("../datas/ontologies_datas")) if 'ontologies_dict' in f]
ontologies_df = pd.concat((pd.read_csv(os.path.join(os.path.abspath("../datas/ontologies_datas/"), f), sep=";", encoding="UTF-8") for f in files), ignore_index=True)

def ontologies_completion(spectrum_list):
    """
        The `ontologies_completion` function enriches a DataFrame containing spectral data with ontological information.

        Parameters:
        - spectrum_list (pd.DataFrame): The DataFrame containing spectral data with 'INCHIKEY' and initial ontology columns.
        - ontologies_df (pd.DataFrame): The DataFrame containing the ontological information to merge.

        Returns:
        - pd.DataFrame: The enriched DataFrame with completed ontology information.

        The function performs the following steps:
        - Initializes columns 'CLASSYFIRE_SUPERCLASS', 'CLASSYFIRE_CLASS', 'CLASSYFIRE_SUBCLASS', 'NPCLASS_PATHWAY',
          'NPCLASS_SUPERCLASS', and 'NPCLASS_CLASS' in `spectrum_list` with default value 'UNKNOWN'.
        - Counts the unique 'INCHIKEY' values in `spectrum_list`.
        - Uses tqdm to create a progress bar for monitoring the process.
        - Merges `spectrum_list` with `ontologies_df` on the 'INCHIKEY' column using a left join.
        - Updates the values in the initial ontology columns with the corresponding values from the merged DataFrame.
        - Removes temporary columns resulted from the merge.
        - Returns the updated DataFrame with completed ontology information.
    """
    # Ajouter les colonnes 'CLASSYFIRE_CLASS' et 'NPCLASSIF_PATHWAY' avec des valeurs par défaut 'UNKNOWN'
    spectrum_list['CLASSYFIRE_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_CLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_SUBCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_PATHWAY'] = "UNKNOWN"
    spectrum_list['NPCLASS_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_CLASS'] = "UNKNOWN"

    # Compter le nombre de INCHIKEY uniques dans spectrum_list
    num_keys = spectrum_list['INCHIKEY'].nunique()

    # Initialise la barre de progression tqdm
    pbar = tqdm(total=num_keys, colour="green", unit=" key", desc="{:>70}".format("updating ontologies"))

    # Fusionner spectrum_list avec ontologies_df sur la colonne 'INCHIKEY'
    completed_df = pd.merge(
        spectrum_list,
        ontologies_df[
            ["INCHIKEY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY",
             "NPCLASS_SUPERCLASS", "NPCLASS_CLASS"]],
        on='INCHIKEY',
        how='left'
    )

    # Mettre à jour les barres de progression
    pbar.update(num_keys)
    pbar.close()

    # Remplacer les valeurs initiales par les valeurs fusionnées
    completed_df['CLASSYFIRE_SUPERCLASS'] = completed_df['CLASSYFIRE_SUPERCLASS_y'].combine_first(
        completed_df['CLASSYFIRE_SUPERCLASS_x'])
    completed_df['CLASSYFIRE_CLASS'] = completed_df['CLASSYFIRE_CLASS_y'].combine_first(
        completed_df['CLASSYFIRE_CLASS_x'])
    completed_df['CLASSYFIRE_SUBCLASS'] = completed_df['CLASSYFIRE_SUBCLASS_y'].combine_first(
        completed_df['CLASSYFIRE_SUBCLASS_x'])
    completed_df['NPCLASS_PATHWAY'] = completed_df['NPCLASS_PATHWAY_y'].combine_first(completed_df['NPCLASS_PATHWAY_x'])
    completed_df['NPCLASS_SUPERCLASS'] = completed_df['NPCLASS_SUPERCLASS_y'].combine_first(
        completed_df['NPCLASS_SUPERCLASS_x'])
    completed_df['NPCLASS_CLASS'] = completed_df['NPCLASS_CLASS_y'].combine_first(completed_df['NPCLASS_CLASS_x'])

    # Supprimer les colonnes temporaires
    completed_df.drop(
        columns=[
            'CLASSYFIRE_SUPERCLASS_x', 'CLASSYFIRE_SUPERCLASS_y',
            'CLASSYFIRE_CLASS_x', 'CLASSYFIRE_CLASS_y',
            'CLASSYFIRE_SUBCLASS_x', 'CLASSYFIRE_SUBCLASS_y',
            'NPCLASS_PATHWAY_x', 'NPCLASS_PATHWAY_y',
            'NPCLASS_SUPERCLASS_x', 'NPCLASS_SUPERCLASS_y',
            'NPCLASS_CLASS_x', 'NPCLASS_CLASS_y'
        ],
        inplace=True
    )

    return completed_df