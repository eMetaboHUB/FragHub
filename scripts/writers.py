from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_json_converted(json_object, original_format):
    """
    Write the given JSON object to a file in the specified format.

    :param json_object: The JSON object to be written.
    :param original_format: The original format of the JSON object.
    :return: None

    Example usage:
        json_object = [...]  # Replace with your JSON object
        original_format = "example"  # Replace with the original format
        write_json_converted(json_object, original_format)

    The method writes the JSON object to a file named "{original_format}_converted.json" in the "../INPUT/JSON/" directory.
    The file is opened in write mode ("w") and encoded using UTF-8. The `ensure_ascii` parameter of the JSON dump function is set to False to preserve non-ASCII characters.

    Note that the method will not write anything if the JSON object is an empty list ([]).
    """
    if json_object != []:
        with open(os.path.join(f"../INPUT/JSON/{original_format}_converted.json"), "w", encoding="UTF-8") as buffer:
            json.dump(json_object, buffer, ensure_ascii=False)

def write_msp(spectrum_list, filename, mode, update):
    """
    Write the given spectrum list to a file in MSP format.

    :param spectrum_list: A list of spectrum strings.
    :param filename: The name of the output file.
    :param mode: The mode to write the file in.
    :return: None
    """
    output_file_path = os.path.join(f"../OUTPUT/MSP/{mode}", filename)

    with tqdm(total=len(spectrum_list), unit=" row", colour="green", desc="{:>80}".format(f"writting {filename}")) as pbar:
        if not update:
            with open(output_file_path, 'w') as f:
                for spectrum in spectrum_list:
                    try:
                        f.write(spectrum)

                        f.write("\n\n")
                    except:
                        continue

                    pbar.update()
        else:
            with open(output_file_path, 'a') as f:
                for spectrum in spectrum_list:
                    try:
                        f.write(spectrum)

                        f.write("\n\n")
                    except:
                        continue

                    pbar.update()



def writting_msp(POS_LC,POS_LC_insilico,POS_GC,POS_GC_insilico,NEG_LC,NEG_LC_insilico,NEG_GC,NEG_GC_insilico, update=False):
    """
    Writes the content of the given parameters to separate MSP files.

    :param POS_LC: The content to write to the "POS_LC.msp" file.
    :param POS_LC_insilico: The content to write to the "POS_LC_insilico.msp" file.
    :param POS_GC: The content to write to the "POS_GC.msp" file.
    :param POS_GC_insilico: The content to write to the "POS_GC_insilico.msp" file.
    :param NEG_LC: The content to write to the "NEG_LC.msp" file.
    :param NEG_LC_insilico: The content to write to the "NEG_LC_insilico.msp" file.
    :param NEG_GC: The content to write to the "NEG_GC.msp" file.
    :param NEG_GC_insilico: The content to write to the "NEG_GC_insilico.msp" file.
    :return: None
    """
    time.sleep(0.1)
    write_msp(POS_LC,"POS_LC.msp", "POS", update)
    del POS_LC
    time.sleep(0.1)
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS", update)
    del POS_LC_insilico
    time.sleep(0.1)
    write_msp(POS_GC, "POS_GC.msp", "POS", update)
    del POS_GC
    time.sleep(0.1)
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS", update)
    del POS_GC_insilico
    time.sleep(0.1)
    write_msp(NEG_LC, "NEG_LC.msp", "NEG", update)
    del NEG_LC
    time.sleep(0.1)
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG", update)
    del NEG_LC_insilico
    time.sleep(0.1)
    write_msp(NEG_GC, "NEG_GC.msp", "NEG", update)
    del NEG_GC
    time.sleep(0.1)
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG", update)
    del NEG_GC_insilico

def write_csv(df, filename, mode, update, first_run):
    """
    :param df: pandas DataFrame object containing the data to be written to CSV.
    :param filename: string representing the name of the output file. The extension ".msp" in the filename will be replaced by ".csv".
    :return: None

    This method writes a pandas DataFrame object to a CSV file. The output file is saved in the "../OUTPUT/CSV/POS" directory with the same name as the input file, but with the extension
    * changed to ".csv". The data is written in chunks of 5000 rows to improve efficiency. The progress of writing is displayed with a progress bar.
    """
    output_file_path = os.path.join(f"../OUTPUT/CSV/{mode}",filename)

    chunk_size = 5000  # Taille de chaque fraction
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))  # Calculer le nombre de fractions

    with tqdm(total=num_chunks, unit=" row", colour="green", desc="{:>80}".format(f"writting {filename}")) as pbar:
        for start in range(0, df.shape[0], chunk_size):
            df_slice = df[start:start + chunk_size]
            if start == 0 and first_run:
                # Écrire les en-têtes pour la première fraction
                df_slice.to_csv(output_file_path, mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)
            else:
                # Append dans le fichier sans écrire les en-têtes pour les autres fractions
                df_slice.to_csv(output_file_path, mode='a', header=False, index=False, sep=";", quotechar='"', encoding="UTF-8")

            # Mettre à jour la barre de progression
            pbar.update()

def writting_csv(POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df,POS_LC_df_insilico,POS_GC_df_insilico,NEG_LC_df_insilico,NEG_GC_df_insilico, first_run, update=False):
    """
    Writes the given dataframes to CSV files with specific file names.

    :param POS_LC_df: DataFrame containing positive LC data
    :param POS_GC_df: DataFrame containing positive GC data
    :param NEG_LC_df: DataFrame containing negative LC data
    :param NEG_GC_df: DataFrame containing negative GC data
    :param POS_LC_df_insilico: DataFrame containing positive LC In Silico data
    :param POS_GC_df_insilico: DataFrame containing positive GC In Silico data
    :param NEG_LC_df_insilico: DataFrame containing negative LC In Silico data
    :param NEG_GC_df_insilico: DataFrame containing negative GC In Silico data
    :return: None
    """
    time.sleep(0.1)
    write_csv(POS_LC_df,"POS_LC.csv","POS", update, first_run)
    del POS_LC_df
    time.sleep(0.1)
    write_csv(POS_GC_df, "POS_GC.csv","POS", update, first_run)
    del POS_GC_df
    time.sleep(0.1)
    write_csv(NEG_LC_df, "NEG_LC.csv","NEG", update, first_run)
    del NEG_LC_df
    time.sleep(0.1)
    write_csv(NEG_GC_df, "NEG_GC.csv","NEG", update, first_run)
    del NEG_GC_df
    time.sleep(0.1)
    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico.csv","POS", update, first_run)
    del POS_LC_df_insilico
    time.sleep(0.1)
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico.csv","POS", update, first_run)
    del POS_GC_df_insilico
    time.sleep(0.1)
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico.csv","NEG", update, first_run)
    del NEG_LC_df_insilico
    time.sleep(0.1)
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico.csv","NEG", update, first_run)
    del NEG_GC_df_insilico