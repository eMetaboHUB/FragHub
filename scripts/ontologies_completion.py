import pandas as pd
import os

global ontologies_df
files = [f for f in os.listdir(os.path.abspath("../datas/ontologies_datas")) if 'ontologies_dict' in f]
ontologies_df = pd.concat((pd.read_csv(os.path.join(os.path.abspath("../datas/ontologies_datas/"), f), sep=";", encoding="UTF-8") for f in files), ignore_index=True)

def ontologies_completion(spectrum_list, progress_callback=None, total_items_callback=None, prefix_callback=None,
                          item_type_callback=None):
    """
    The `ontologies_completion` function enriches a DataFrame containing spectral data with ontological information.

    Parameters:
    - spectrum_list (pd.DataFrame): The DataFrame containing spectral data with 'INCHIKEY' and initial ontology columns.
    - progress_callback (callable, optional): Function to update progress during processing.
    - total_items_callback (callable, optional): Function to set the total number of items to process.
    - prefix_callback (callable, optional): Function to describe the task being performed.
    - item_type_callback (callable, optional): Function to define the item type being processed.

    Returns:
    - pd.DataFrame: The enriched DataFrame with completed ontology information.
    """
    # Ajouter les colonnes avec des valeurs par défaut 'UNKNOWN'
    spectrum_list['CLASSYFIRE_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_CLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_SUBCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_PATHWAY'] = "UNKNOWN"
    spectrum_list['NPCLASS_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_CLASS'] = "UNKNOWN"

    # Compter le nombre de INCHIKEY uniques dans spectrum_list
    num_keys = spectrum_list['INCHIKEY'].nunique()

    # Initialiser le suivi avec les callbacks
    if prefix_callback:
        prefix_callback("updating ontologies:")

    if item_type_callback:
        item_type_callback("rows")

    if total_items_callback:
        total_items_callback(num_keys, 0)  # Définir le total au départ

    # Fusionner spectrum_list avec ontologies_df sur 'INCHIKEY'
    completed_df = pd.merge(
        spectrum_list,
        ontologies_df[
            ["INCHIKEY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY",
             "NPCLASS_SUPERCLASS", "NPCLASS_CLASS"]
        ],
        on='INCHIKEY',
        how='left'
    )

    # Simuler mise à jour clé par clé pour la progression
    processed_keys = 0
    for _ in range(num_keys):
        processed_keys += 1
        if progress_callback:
            progress_callback(processed_keys)

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
