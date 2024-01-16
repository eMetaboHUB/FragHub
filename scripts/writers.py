from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_msp(spectrum_list, filename, mode, update):
    """
    Write the given spectrum list to a file in MSP format.

    :param spectrum_list: A list of spectrum strings.
    :param filename: The name of the output file.
    :param mode: The mode to write the file in.
    :return: None
    """
    output_file_path = os.path.join(f"../OUTPUT/MSP/{mode}", filename)

    with tqdm(total=len(spectrum_list), unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:
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

    with tqdm(total=num_chunks, unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:
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

def write_json(df, filename, mode):
    """
    :param df: pandas DataFrame object containing the data to be written to JSON.
    :param filename: string representing the name of the output file. The extension ".csv" in the filename will be replaced by ".json".
    :return: None
    This method writes a pandas DataFrame object to a JSON array in file. The output file is saved in the "../OUTPUT/JSON/{mode}" directory with the same name as the input file,
    but with the extension changed to ".json".
    """
    output_file_path = os.path.join(f"../OUTPUT/JSON/{mode}", filename)

    # Convert the DataFrame to a list of dict records
    json_records = df.to_dict('records')

    with open(output_file_path, 'w') as f:
        pbar = tqdm(total=len(json_records), unit=" row", colour="green", desc="{:>70}".format(f"writing {filename}"))
        f.write('[\n')
        for i, record in enumerate(json_records):
            f.write('\t')
            json.dump(record, f)
            if i < len(json_records) - 1:  # Avoids adding a comma after the last record
                f.write(',\n')
            pbar.update()
        f.write('\n]')

def writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico):
    """
    Write JSON files for the given data frames.

    :param POS_LC_df: Positive LC data frame
    :param POS_GC_df: Positive GC data frame
    :param NEG_LC_df: Negative LC data frame
    :param NEG_GC_df: Negative GC data frame
    :param POS_LC_df_insilico: Positive LC In Silico data frame
    :param POS_GC_df_insilico: Positive GC In Silico data frame
    :param NEG_LC_df_insilico: Negative LC In Silico data frame
    :param NEG_GC_df_insilico: Negative GC In Silico data frame
    :return: None
    """
    time.sleep(0.1)
    write_json(POS_LC_df, "POS_LC.json", "POS")
    del POS_LC_df
    time.sleep(0.1)
    write_json(POS_GC_df, "POS_GC.json", "POS")
    del POS_GC_df
    time.sleep(0.1)
    write_json(NEG_LC_df, "NEG_LC.json", "NEG")
    del NEG_LC_df
    time.sleep(0.1)
    write_json(NEG_GC_df, "NEG_GC.json", "NEG")
    del NEG_GC_df
    time.sleep(0.1)
    write_json(POS_LC_df_insilico, "POS_LC_In_Silico.json", "POS")
    del POS_LC_df_insilico
    time.sleep(0.1)
    write_json(POS_GC_df_insilico, "POS_GC_In_Silico.json", "POS")
    del POS_GC_df_insilico
    time.sleep(0.1)
    write_json(NEG_LC_df_insilico, "NEG_LC_In_Silico.json", "NEG")
    del NEG_LC_df_insilico
    time.sleep(0.1)
    write_json(NEG_GC_df_insilico, "NEG_GC_In_Silico.json", "NEG")
    del NEG_GC_df_insilico
