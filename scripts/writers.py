from tqdm import tqdm
import pandas as pd
import numpy as np
import time
import os
import re

def write_msp(spectrum_list, filename, mode):
    """
    Write the given spectrum list to a file in MSP format.

    :param spectrum_list: A list of spectrum strings.
    :param filename: The name of the output file.
    :param mode: The mode to write the file in.
    :return: None
    """
    time.sleep(0.1)
    print(f"-- {filename.replace('.msp', '.csv')} --")

    output_file_path = os.path.join(f"../OUTPUT/MSP/{mode}", filename)

    with tqdm(total=len(spectrum_list), unit=" row", colour="green", desc="{:>25}".format("writting")) as pbar:
        with open(output_file_path, 'w') as f:
            for spectrum in spectrum_list:
                try:
                    f.write(spectrum)

                    f.write("\n\n\n")
                except:
                    continue

                pbar.update()


def writting_msp(POS_LC,POS_LC_insilico,POS_GC,POS_GC_insilico,NEG_LC,NEG_LC_insilico,NEG_GC,NEG_GC_insilico):
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
    write_msp(POS_LC,"POS_LC.msp", "POS")
    del POS_LC
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS")
    del POS_LC_insilico
    write_msp(POS_GC, "POS_GC.msp", "POS")
    del POS_GC
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS")
    del POS_GC_insilico
    write_msp(NEG_LC, "NEG_LC.msp", "NEG")
    del NEG_LC
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG")
    del NEG_LC_insilico
    write_msp(NEG_GC, "NEG_GC.msp", "NEG")
    del NEG_GC
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG")
    del NEG_GC_insilico

def write_csv(df, filename, mode):
    """
    :param df: pandas DataFrame object containing the data to be written to CSV.
    :param filename: string representing the name of the output file. The extension ".msp" in the filename will be replaced by ".csv".
    :return: None

    This method writes a pandas DataFrame object to a CSV file. The output file is saved in the "../OUTPUT/CSV/POS" directory with the same name as the input file, but with the extension
    * changed to ".csv". The data is written in chunks of 5000 rows to improve efficiency. The progress of writing is displayed with a progress bar.
    """
    time.sleep(0.1)
    print(f"-- {filename.replace('.msp','.csv')} --")

    output_file_path = os.path.join(f"../OUTPUT/CSV/{mode}",filename)

    chunk_size = 5000  # Taille de chaque fraction
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))  # Calculer le nombre de fractions

    with tqdm(total=num_chunks, unit=" row", colour="green", desc="{:>25}".format("writting")) as pbar:
        for start in range(0, df.shape[0], chunk_size):
            df_slice = df[start:start + chunk_size]
            if start == 0:
                # Écrire les en-têtes pour la première fraction
                df_slice.to_csv(output_file_path, mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)
            else:
                # Append dans le fichier sans écrire les en-têtes pour les autres fractions
                df_slice.to_csv(output_file_path, mode='a', header=False, index=False, sep=";", quotechar='"', encoding="UTF-8")

            # Mettre à jour la barre de progression
            pbar.update()

def writting_csv(POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df,POS_LC_df_insilico,POS_GC_df_insilico,NEG_LC_df_insilico,NEG_GC_df_insilico):
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
    write_csv(POS_LC_df,"POS_LC_clean.csv","POS")
    del POS_LC_df
    write_csv(POS_GC_df, "POS_GC_clean.csv","POS")
    del POS_GC_df
    write_csv(NEG_LC_df, "NEG_LC_clean.csv","NEG")
    del NEG_LC_df
    write_csv(NEG_GC_df, "NEG_GC_clean.csv","NEG")
    del NEG_GC_df

    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico_clean.csv","POS")
    del POS_LC_df_insilico
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico_clean.csv","POS")
    del POS_GC_df_insilico
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico_clean.csv","NEG")
    del NEG_LC_df_insilico
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico_clean.csv","NEG")
    del NEG_GC_df_insilico