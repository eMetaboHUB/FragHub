import pandas as pd
import ast

import time


def remove_duplicatas(spectrum_list, progress_callback=None, total_items_callback=None, prefix_callback=None,
                      item_type_callback=None):
    """
    Remove duplicate entries from a given spectrum list based on the maximum
    'row_size' for each unique SPLASH identifier. The function computes a temporary
    'row_size' for each row, removes duplicates by keeping the entry with the largest
    'row_size' for each SPLASH, and returns a list of dictionaries representing the cleaned data.

    Parameters:
    spectrum_list : DataFrame
        A pandas DataFrame where each row represents a spectrum, and each spectrum
        must have a 'SPLASH' identifier used to identify duplicates in the DataFrame.
    progress_callback : callable, optional
        A function to report progress (processed items).
    total_items_callback : callable, optional
        A function to report the total number of items to process.
    prefix_callback : callable, optional
        A function to set the prefix for the operation.
    item_type_callback : callable, optional
        A function to specify the type of items.

    Returns:
    list of dict
        A list of dictionaries with duplicate spectrums removed, retaining only the
        spectrum with the largest 'row_size' for each unique SPLASH identifier.
    """
    total_items = len(spectrum_list)

    # Définir un préfixe via le callback, si fourni
    if prefix_callback:
        prefix_callback("removing duplicates")

    # Spécifier le type d'éléments via le callback
    if item_type_callback:
        item_type_callback("spectra")

    # Définir le total via le callback, si fourni
    if total_items_callback:
        total_items_callback(len(spectrum_list), 0)  # total = nombre initial de spectres

    # Calculer la taille des lignes (en caractères) en ajoutant une nouvelle colonne temporaire 'row_size'
    spectrum_list['row_size'] = spectrum_list.apply(lambda row: row.astype(str).map(len).sum(), axis=1)

    # Supprimer les doublons en gardant la ligne avec la plus grande 'row_size' pour chaque SPLASH
    unique_splash_indices = spectrum_list.groupby('SPLASH')['row_size'].idxmax()
    spectrum_list = spectrum_list.loc[unique_splash_indices]

    # Supprimer la colonne temporaire 'row_size'
    spectrum_list = spectrum_list.drop(columns=['row_size'])

    # Mettre à jour la progression après le traitement (et passer à 100 % à la fin)
    if progress_callback:
        progress_callback(total_items)  # Mise à jour avec le nombre final de spectres
        progress_callback(100)  # Forcer la progression à 100 %

    # Convertir en liste de dictionnaires pour la sortie
    spectrum_list = spectrum_list.to_dict(orient='records')

    return spectrum_list
