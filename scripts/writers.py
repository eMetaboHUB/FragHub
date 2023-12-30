from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import re

def writting_msp(clean_msp_path,POS_LC,POS_GC,NEG_LC,NEG_GC,POS_LC_In_Silico,POS_GC_In_Silico,NEG_LC_In_Silico,NEG_GC_In_Silico):
    """
    :param clean_msp_path: The path where the cleaned MSP files will be saved.
    :param POS_LC: A list of strings representing the positive LC data.
    :param POS_GC: A list of strings representing the positive GC data.
    :param NEG_LC: A list of strings representing the negative LC data.
    :param NEG_GC: A list of strings representing the negative GC data.
    :param POS_LC_In_Silico: A list of strings representing the positive LC In Silico data.
    :param POS_GC_In_Silico: A list of strings representing the positive GC In Silico data.
    :param NEG_LC_In_Silico: A list of strings representing the negative LC In Silico data.
    :param NEG_GC_In_Silico: A list of strings representing the negative GC In Silico data.
    :return: None

    This method writes the cleaned MSP files to the specified `clean_msp_path` location. The cleaned data is obtained by joining the strings in the respective input lists and applying regular
    * expression substitutions.

    The cleaned MSP files are saved with the following names and paths:
    - 'POS_LC_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'POS_GC_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'NEG_LC_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'NEG_GC_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'POS_LC_In_Silico_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'POS_GC_In_Silico_clean.msp' in the 'POS' folder under `clean_msp_path`
    - 'NEG_LC_In_Silico_clean.msp' in the 'NEG' folder under `clean_msp_path`
    - 'NEG_GC_In_Silico_clean.msp' in the 'NEG' folder under `clean_msp_path`
    """
    POS_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(POS_LC))
    del POS_LC
    POS_GC_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_GC))
    del POS_GC
    NEG_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_LC))
    del NEG_LC
    NEG_GC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_GC))
    del NEG_GC

    POS_LC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_LC_In_Silico))
    del POS_LC_In_Silico
    POS_GC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_GC_In_Silico))
    del POS_GC_In_Silico
    NEG_LC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(NEG_LC_In_Silico))
    del NEG_LC_In_Silico
    NEG_GC_In_Silico_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(NEG_GC_In_Silico))
    del NEG_GC_In_Silico

    with open(os.path.join(clean_msp_path,"POS/POS_LC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_LC_FULL)
    del POS_LC_FULL
    with open(os.path.join(clean_msp_path,"POS/POS_GC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_GC_FULL)
    del POS_GC_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_LC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_LC_FULL)
    del NEG_LC_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_GC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_GC_FULL)
    del NEG_GC_FULL

    with open(os.path.join(clean_msp_path,"POS/POS_LC_In_Silico_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_LC_In_Silico_FULL)
    del POS_LC_In_Silico_FULL
    with open(os.path.join(clean_msp_path,"POS/POS_GC_In_Silico_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_GC_In_Silico_FULL)
    del POS_GC_In_Silico_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_LC_In_Silico_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_LC_In_Silico_FULL)
    del NEG_LC_In_Silico_FULL
    with open(os.path.join(clean_msp_path, "NEG/NEG_GC_In_Silico_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_GC_In_Silico_FULL)
    del NEG_GC_In_Silico_FULL

def write_csv(df,filename):
    """
    :param df: pandas DataFrame object containing the data to be written to CSV.
    :param filename: string representing the name of the output file. The extension ".msp" in the filename will be replaced by ".csv".
    :return: None

    This method writes a pandas DataFrame object to a CSV file. The output file is saved in the "../OUTPUT/CSV/POS" directory with the same name as the input file, but with the extension
    * changed to ".csv". The data is written in chunks of 5000 rows to improve efficiency. The progress of writing is displayed with a progress bar.
    """
    print(f"-- {filename.replace('.msp','.csv')} --")

    output_file_path = os.path.join("../OUTPUT/CSV/POS",filename)

    chunk_size = 5000  # Taille de chaque fraction
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))  # Calculer le nombre de fractions

    with tqdm(total=num_chunks, unit=" row", colour="green", desc="\t    writting") as pbar:
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
    write_csv(POS_LC_df,"POS_LC_clean.csv")
    del POS_LC_df
    write_csv(POS_GC_df, "POS_GC_clean.csv")
    del POS_GC_df
    write_csv(NEG_LC_df, "NEG_LC_clean.csv")
    del NEG_LC_df
    write_csv(NEG_GC_df, "NEG_GC_clean.csv")
    del NEG_GC_df

    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico_clean.csv")
    del POS_LC_df_insilico
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico_clean.csv")
    del POS_GC_df_insilico
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico_clean.csv")
    del NEG_LC_df_insilico
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico_clean.csv")
    del NEG_GC_df_insilico