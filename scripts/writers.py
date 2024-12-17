from set_projects import parameters_dict
from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_msp(spectrum_list, filename, mode, update, profile_name, progress_callback=None, total_items_callback=None,
              prefix_callback=None, item_type_callback=None):
    """
    Write MSP file.
    :param spectrum_list: A list of spectra to be written to the MSP file.
    :param filename: The name of the output file.
    :param mode: The mode of writing (either 'w' for write or 'a' for append).
    :param update: Boolean value indicating whether to update the existing file or create a new one.
    :param profile_name: The name of the profile being used.
    :param progress_callback: Callback for reporting progress.
    :param total_items_callback: Callback for setting the total number of items.
    :param prefix_callback: Callback to set the prefix for context/task.
    :param item_type_callback: Callback for specifying the type of items (e.g., "row").
    :return: None
    """

    # Gérer le chemin de sortie
    output_file_path = os.path.join(parameters_dict["output_directory"], f"{profile_name}/MSP/{mode}", filename)

    # Donne le contexte du traitement avec le prefix_callback
    if prefix_callback:
        prefix_callback(f"Writing {filename} to MSP")

    # Indiquer le type d'éléments avec le item_type_callback
    if item_type_callback:
        item_type_callback("spectra")

    # Initialiser le total d'éléments avec le total_items_callback
    total_spectra = len(spectrum_list)
    if total_items_callback:
        total_items_callback(total_spectra, 0)  # Initialisation à 0 progressé

    # Ouvrir le fichier pour écrire les spectres
    file_mode = 'a' if update else 'w'
    with open(output_file_path, file_mode, encoding="UTF-8") as f:
        # Parcourir chaque spectre dans la liste
        for index, spectrum in enumerate(spectrum_list):
            try:
                # Écrire le spectre dans le fichier
                f.write(spectrum)
                f.write("\n\n")

                # Mettre à jour les callbacks après chaque écriture
                if progress_callback:
                    progress_callback(index + 1)  # Signal progress (ligne actuelle)
            except Exception as e:
                # En cas d'erreur, ignorer le spectre, mais signaler les problèmes si nécessaire
                print(f"Error writing spectrum {index}: {e}")
                continue


def writting_msp(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico, profile_name, update=False, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    This function writes MSP files for positive and negative LC and GC data. It takes in various inputs for different types of data, a profile_name which is the name of the profile and an optional update flag which defaults to False when not provided.

    Args:
      POS_LC: This parameter represents the positive LC data.
      POS_LC_insilico: This parameter represents the insilico data for positive LC.
      POS_GC: This parameter represents the positive GC data.
      POS_GC_insilico: This parameter represents the insilico data for positive GC.
      NEG_LC: This parameter represents the negative LC data.
      NEG_LC_insilico: This parameter represents the insilico data for negative LC.
      NEG_GC: This parameter represents the negative GC data.
      NEG_GC_insilico: This parameter represents the insilico data for negative GC.
      profile_name: This is the name of the profile.
      update: This is an optional parameter which is set to False by default. If set to True, it results in an update operation.

    Returns: None
    """

    # sleep for 0.1 second to avoid too many computations
    time.sleep(0.1)
    # write the positive LC data to the "POS_LC.msp" file
    write_msp(POS_LC, "POS_LC.msp", "POS", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # delete the positive LC data to free up memory
    del POS_LC
    # repeat the process for other types of data
    time.sleep(0.1)
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_insilico
    time.sleep(0.1)
    write_msp(POS_GC, "POS_GC.msp", "POS", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC
    time.sleep(0.1)
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_insilico
    time.sleep(0.1)
    write_msp(NEG_LC, "NEG_LC.msp", "NEG", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC
    time.sleep(0.1)
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_insilico
    time.sleep(0.1)
    write_msp(NEG_GC, "NEG_GC.msp", "NEG", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC
    time.sleep(0.1)
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG", update, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_insilico


def write_csv(df, filename, mode, update, first_run, profile_name, progress_callback=None, total_items_callback=None,
              prefix_callback=None, item_type_callback=None):
    """
    This function writes a pandas DataFrame to a CSV file in chunks.
    :param df: DataFrame - the data to be written to the CSV file.
    :param filename: str - the name of the CSV file to write the data to.
    :param mode: str - the mode in which the file should be opened ('w' - write, 'a' - append).
    :param update: Bool - a flag to indicate whether the file is updated or created from scratch.
    :param first_run: Bool - indicates whether it's the first run (write the headers if True).
    :param profile_name: str - a name related to the CSV file for output organization.
    :param progress_callback: Callable, function to track progress of row writing.
    :param total_items_callback: Callable, function to set the total number of rows.
    :param prefix_callback: Callable, function to indicate the current task's description context.
    :param item_type_callback: Callable, function to indicate the type of elements being processed.
    :return: None.
    """
    # Replace newline characters with semicolons in the PEAKS_LIST column
    if 'PEAKS_LIST' in df.columns:
        df['PEAKS_LIST'] = df['PEAKS_LIST'].str.replace('\n', ';', regex=False)

    # Construct the file path dynamically
    output_file_path = os.path.join(parameters_dict["output_directory"], f"{profile_name}/CSV/{mode}", filename)

    # Define chunk size for writing DataFrame in parts
    chunk_size = 5000

    # Calculate the total number of chunks
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))

    # Set up task description
    if prefix_callback:
        prefix_callback(f"Writing {filename} to CSV")

    # Indicate item type
    if item_type_callback:
        item_type_callback("DataFrame Rows")

    # Define total items expected
    total_rows = len(df)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Initialize total and completed count

    # Process the DataFrame in chunks and write each chunk to the file
    for chunk_index, start in enumerate(range(0, total_rows, chunk_size)):
        # Select the slice of the DataFrame for the current chunk
        df_slice = df.iloc[start:start + chunk_size]

        # Write headers only on the first chunk and during the first run
        if start == 0 and first_run:
            df_slice.to_csv(output_file_path, mode='w', sep="\t", quotechar='"', encoding="UTF-8", index=False)
        else:
            # Append the chunk without writing headers
            df_slice.to_csv(output_file_path, mode='a', sep="\t", quotechar='"', encoding="UTF-8", index=False,
                            header=False)

        # Notify progress of written chunks (via progress_callback)
        if progress_callback:
            processed_rows = min((chunk_index + 1) * chunk_size, total_rows)
            progress_callback(processed_rows)  # Update with processed row count


def writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, profile_name, update=False, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write data to CSV files.

    :param POS_LC_df: Pandas DataFrame containing Positive LC data
    :param POS_GC_df: Pandas DataFrame containing Positive GC data
    :param NEG_LC_df: Pandas DataFrame containing Negative LC data
    :param NEG_GC_df: Pandas DataFrame containing Negative GC data
    :param POS_LC_df_insilico: Pandas DataFrame containing Positive LC data (In Silico)
    :param POS_GC_df_insilico: Pandas DataFrame containing Positive GC data (In Silico)
    :param NEG_LC_df_insilico: Pandas DataFrame containing Negative LC data (In Silico)
    :param NEG_GC_df_insilico: Pandas DataFrame containing Negative GC data (In Silico)
    :param first_run: Boolean indicating if it's the first run
    :param profile_name: Name of the profile
    :param update: Boolean indicating if it's an update
    :return: None
    """

    # Pause execution for 0.1 seconds
    time.sleep(0.1)
    # Call write_csv for Positive LC data and write it to POS_LC.csv file
    write_csv(POS_LC_df, "POS_LC.csv", "POS", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Clear memory used by POS_LC_df dataframe
    del POS_LC_df

    # Similar comments can be written for the rest of the calls and deletions below
    # Just replace "POS_LC_df" and "POS_LC.csv" with appropriate values

    time.sleep(0.1)
    write_csv(POS_GC_df, "POS_GC.csv", "POS", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df

    time.sleep(0.1)
    write_csv(NEG_LC_df, "NEG_LC.csv", "NEG", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df

    time.sleep(0.1)
    write_csv(NEG_GC_df, "NEG_GC.csv", "NEG", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df

    time.sleep(0.1)
    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico.csv", "POS", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_df_insilico

    time.sleep(0.1)
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico.csv", "POS", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico.csv", "NEG", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df_insilico

    time.sleep(0.1)
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico.csv", "NEG", update, first_run, profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df_insilico


def write_json(df, filename, mode, profile_name, progress_callback=None, total_items_callback=None,
               prefix_callback=None, item_type_callback=None):
    """
    Write DataFrame to a JSON file.
    :param df: The DataFrame to be written.
    :type df: pandas.DataFrame
    :param filename: The name of the output file.
    :type filename: str
    :param mode: The mode for opening the output file. e.g., 'w', 'a'.
    :type mode: str
    :param profile_name: The name of the profile.
    :type profile_name: str
    :param progress_callback: Callback to report progress.
    :param total_items_callback: Callback to set total items count.
    :param prefix_callback: Callback to define a prefix or task description.
    :param item_type_callback: Callback to specify the type of items being processed.
    :return: None
    """

    # Définir le chemin du fichier de sortie
    output_file_path = os.path.join(parameters_dict["output_directory"], f"{profile_name}/JSON/{mode}", filename)

    # Convertir le DataFrame en une liste de dictionnaires
    json_records = df.to_dict('records')

    # Si un prefix_callback est défini, indiquer la tâche actuelle
    if prefix_callback:
        prefix_callback(f"Writing {filename} to JSON")

    # Définir le type d'éléments avec le callback item_type_callback
    if item_type_callback:
        item_type_callback("rows")

    # Déclarer le nombre total d'enregistrements à traiter
    total_records = len(json_records)
    if total_items_callback:
        total_items_callback(total_records, 0)  # Initialiser les éléments traités à 0

    # Ouvrir le fichier pour écrire
    with open(output_file_path, 'w') as f:
        # Écrire le crochet ouvrant pour démarrer le tableau JSON
        f.write('[\n')

        # Parcourir les enregistrements convertis
        for i, record in enumerate(json_records):
            # Écrire une tabulation pour l'indentation
            f.write('\t')

            # Écrire le dictionnaire sous forme de chaîne JSON
            json.dump(record, f)

            # Ajouter une virgule à la fin, sauf pour le dernier élément
            if i < total_records - 1:
                f.write(',\n')

            # S'il y a un callback de progression, l'appeler pour mettre à jour la progression
            if progress_callback:
                progress_callback(i + 1)  # Mise à jour avec l'index actuel (commençant à 1)

        # Écrire le crochet fermant pour terminer le tableau JSON
        f.write('\n]')


def writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, profile_name, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Write dataframes to JSON files.
    :param POS_LC_df: Positive LC dataframe
    :param POS_GC_df: Positive GC dataframe
    :param NEG_LC_df: Negative LC dataframe
    :param NEG_GC_df: Negative GC dataframe
    :param POS_LC_df_insilico: Positive LC In Silico dataframe
    :param POS_GC_df_insilico: Positive GC In Silico dataframe
    :param NEG_LC_df_insilico: Negative LC In Silico dataframe
    :param NEG_GC_df_insilico: Negative GC In Silico dataframe
    :param profile_name: Name of the profile
    :return: None
    """
    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive LC DataFrame to a JSON file named "POS_LC.json"
    write_json(POS_LC_df, "POS_LC.json", "POS", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Deleting the DataFrame from memory as it's no longer needed
    del POS_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive GC DataFrame to a JSON file named "POS_GC.json"
    write_json(POS_GC_df, "POS_GC.json", "POS", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative LC DataFrame to a JSON file named "NEG_LC.json"
    write_json(NEG_LC_df, "NEG_LC.json", "NEG", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative GC DataFrame to a JSON file named "NEG_GC.json"
    write_json(NEG_GC_df, "NEG_GC.json", "NEG", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_GC_df

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico LC DataFrame to a JSON file named "POS_LC_In_Silico.json"
    write_json(POS_LC_df_insilico, "POS_LC_In_Silico.json", "POS", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Positive In Silico GC DataFrame to a JSON file named "POS_GC_In_Silico.json"
    write_json(POS_GC_df_insilico, "POS_GC_In_Silico.json", "POS", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del POS_GC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico LC DataFrame to a JSON file named "NEG_LC_In_Silico.json"
    write_json(NEG_LC_df_insilico, "NEG_LC_In_Silico.json", "NEG", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    del NEG_LC_df_insilico

    time.sleep(0.1)  # Providing a short delay to ensure smooth execution of next command

    # Write the Negative In Silico GC DataFrame to a JSON file named "NEG_GC_In_Silico.json"
    write_json(NEG_GC_df_insilico, "NEG_GC_In_Silico.json", "NEG", profile_name, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Deleting the DataFrame from memory as it's no longer needed
    del NEG_GC_df_insilico
