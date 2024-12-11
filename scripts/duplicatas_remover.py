from tqdm import tqdm
import pandas as pd
import ast

import time

def remove_duplicatas(spectrum_list):
    """
    Remove duplicate entries from a given spectrum list based on the maximum
    'row_size' for each unique SPLASH identifier. The input is a DataFrame where
    each row represents a spectrum and contains a unique identifier called SPLASH.
    The function computes a temporary 'row_size' for each row, removes duplicates
    by keeping the entry with the largest 'row_size' for each SPLASH, and returns
    a list of dictionaries representing the cleaned data.

    Parameters:
    spectrum_list : DataFrame
        A pandas DataFrame where each row represents a spectrum, and each spectrum
        must have a 'SPLASH' identifier used to identify duplicates in the DataFrame.

    Returns:
    list of dict
        A list of dictionaries with duplicate spectrums removed, retaining only the
        spectrum with the largest 'row_size' for each unique SPLASH identifier.
    """
    # Calculer la taille en bytes de chaque ligne en ajoutant une nouvelle colonne temporaire 'row_size'
    spectrum_list['row_size'] = spectrum_list.apply(lambda row: row.astype(str).map(len).sum(), axis=1)

    # Configurer tqdm pour suivre le processus de suppression
    desc = "{:>70}".format("removing duplicates")
    with tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc=desc) as pbar:
        # Supprimer les doublons en gardant la ligne où 'row_size' est maximale pour chaque SPLASH
        spectrum_list = spectrum_list.loc[spectrum_list.groupby('SPLASH')['row_size'].idxmax()]
        pbar.update(len(spectrum_list))  # Avancement

    # Supprimer la colonne temporaire 'row_size'
    spectrum_list = spectrum_list.drop(columns=['row_size'])

    spectrum_list = spectrum_list.to_dict(orient='records')

    return spectrum_list
