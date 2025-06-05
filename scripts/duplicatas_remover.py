import scripts.deletion_report
import pandas as pd
import os

def remove_duplicatas(spectrum_list, output_directory, progress_callback=None, total_items_callback=None,
                      prefix_callback=None, item_type_callback=None):
    """
    Remove duplicate entries from a given spectrum list based on the maximum 'row_size' for each unique
    combination of SPLASH and INCHIKEY. Entries with empty INCHIKEY are not considered for deduplication.

    Parameters:
    -----------
    spectrum_list : DataFrame
        A pandas DataFrame where each row represents a spectrum. Each spectrum must have
        'SPLASH' and 'INCHIKEY' columns used to identify duplicates.

    output_directory : str
        The base directory where `DELETED_SPECTRUMS` and the CSV file will be stored.

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
        spectrum with the largest 'row_size' for each combination of SPLASH and INCHIKEY.
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

    # Séparer les entrées avec INCHIKEY vide ou non
    empty_inchikey_df = spectrum_list[spectrum_list['INCHIKEY'].isna() | (spectrum_list['INCHIKEY'].str.strip() == "")]
    non_empty_inchikey_df = spectrum_list[~(spectrum_list['INCHIKEY'].isna() | (spectrum_list['INCHIKEY'].str.strip() == ""))]

    # Identifier les indices uniques pour SPLASH et INCHIKEY (non vides)
    unique_indices = non_empty_inchikey_df.groupby(['SPLASH', 'INCHIKEY'])['row_size'].idxmax()

    # Combiner les indices à garder : uniques et entrées sans INCHIKEY
    all_indices_to_keep = set(unique_indices).union(set(empty_inchikey_df.index))

    # Identifier les doublons à supprimer
    all_indices = set(spectrum_list.index)
    indices_to_delete = list(all_indices - all_indices_to_keep)

    # Extraire les doublons à supprimer
    deleted_spectra = spectrum_list.loc[indices_to_delete].copy()

    # Supprimer la colonne temporaire 'row_size' avant de les ajouter à la liste des doublons
    deleted_spectra = deleted_spectra.drop(columns=['row_size'])
    deleted_spectra['DELETION_REASON'] = "spectrum deleted because it's a duplicate (SPLASH + INCHIKEY)"

    # Créer le répertoire pour stocker les spectres supprimés
    deleted_spectrums_dir = os.path.join(output_directory, 'DELETED_SPECTRUMS')
    os.makedirs(deleted_spectrums_dir, exist_ok=True)

    # Écrire les doublons supprimés dans un fichier CSV
    deleted_spectra_file = os.path.join(deleted_spectrums_dir, 'duplicatas_removed.csv')
    deleted_spectra.to_csv(deleted_spectra_file, sep='\t', index=False, quotechar='"')

    del deleted_spectra

    # Garder uniquement les spectres uniques (conversion du set en liste)
    spectrum_list = spectrum_list.loc[list(all_indices_to_keep)]

    # Supprimer la colonne temporaire 'row_size'
    spectrum_list = spectrum_list.drop(columns=['row_size'])

    # Mettre à jour la progression, si nécessaire
    if progress_callback:
        progress_callback(total_items)
        progress_callback(100)

    # Convertir en liste de dictionnaires
    spectrum_list = spectrum_list.to_dict(orient='records')

    # Mettre à jour le rapport de suppression
    scripts.deletion_report.duplicatas_removed = total_items - len(spectrum_list)

    return spectrum_list
