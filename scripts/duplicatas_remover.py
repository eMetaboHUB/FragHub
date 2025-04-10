import deletion_report
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
    -----------
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
    --------
    list of dict
        A list of dictionaries with duplicate spectrums removed, retaining only the
        spectrum with the largest 'row_size' for each unique SPLASH identifier.
    """

    total_items = len(spectrum_list)

    if prefix_callback:
        prefix_callback("Removing duplicates:")

    if item_type_callback:
        item_type_callback("spectra")

    if total_items_callback:
        total_items_callback(total_items, 0)

    # Calculer la taille des lignes en caractères
    spectrum_list['row_size'] = spectrum_list.apply(lambda row: row.astype(str).map(len).sum(), axis=1)

    # Identifier les SPLASH uniques avec la plus grande 'row_size'
    unique_splash_indices = spectrum_list.groupby('SPLASH')['row_size'].idxmax()

    # Identifier les doublons à supprimer
    all_indices = set(spectrum_list.index)
    indices_to_keep = set(unique_splash_indices)
    indices_to_delete = list(all_indices - indices_to_keep)

    # Extraire les doublons à supprimer
    deleted_spectra = spectrum_list.loc[indices_to_delete].copy()

    # Supprimer la colonne temporaire 'row_size' avant de les ajouter à la liste des doublons
    deleted_spectra = deleted_spectra.drop(columns=['row_size'])
    deleted_spectra['DELETION_REASON'] = "spectrum deleted because it's a duplicatas"

    # Ajouter les spectres supprimés à la liste deletion_report.deleted_spectrum_list
    deletion_report.deleted_spectrum_list.extend(deleted_spectra.to_dict(orient='records'))

    # Garder uniquement les spectres uniques
    spectrum_list = spectrum_list.loc[unique_splash_indices]

    # Supprimer la colonne temporaire 'row_size'
    spectrum_list = spectrum_list.drop(columns=['row_size'])

    # Mettre à jour la progression, si nécessaire
    if progress_callback:
        progress_callback(total_items)
        progress_callback(100)

    # Convertir en liste de dictionnaires
    spectrum_list = spectrum_list.to_dict(orient='records')

    # Mettre à jour le rapport de suppression
    deletion_report.duplicatas_removed = total_items - len(spectrum_list)

    return spectrum_list
