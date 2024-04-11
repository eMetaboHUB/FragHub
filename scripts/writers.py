from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import time
import os
import re

def write_msp(spectrum_list, filename, mode, update, profile_name):
    """
    Write MSP file.

    :param spectrum_list: A list of spectra to be written to the MSP file.
    :param filename: The name of the output file.
    :param mode: The mode of writing (either 'w' for write or 'a' for append).
    :param update: Boolean value indicating whether to update the existing file or create a new one.
    :param profile_name: The name of the profile being used.
    :return: None
    """
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/MSP/{mode}", filename)

    with tqdm(total=len(spectrum_list), unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:
        if not update:
            with open(output_file_path, 'w', encoding="UTF-8") as f:
                for spectrum in spectrum_list:
                    try:
                        f.write(spectrum)

                        f.write("\n\n")
                    except:
                        continue

                    pbar.update()
        else:
            with open(output_file_path, 'a', encoding="UTF-8") as f:
                for spectrum in spectrum_list:
                    try:
                        f.write(spectrum)

                        f.write("\n\n")
                    except:
                        continue

                    pbar.update()



def writting_msp(POS_LC,POS_LC_insilico,POS_GC,POS_GC_insilico,NEG_LC,NEG_LC_insilico,NEG_GC,NEG_GC_insilico, profile_name, update=False):
    """
    Write MSP files for positive and negative LC and GC data.

    :param POS_LC: positive LC data
    :param POS_LC_insilico: positive LC data insilico
    :param POS_GC: positive GC data
    :param POS_GC_insilico: positive GC data insilico
    :param NEG_LC: negative LC data
    :param NEG_LC_insilico: negative LC data insilico
    :param NEG_GC: negative GC data
    :param NEG_GC_insilico: negative GC data insilico
    :param profile_name: name of the profile
    :param update: (optional) update flag, defaults to False
    :return: None
    """
    time.sleep(0.1)
    write_msp(POS_LC,"POS_LC.msp", "POS", update, profile_name)
    del POS_LC
    time.sleep(0.1)
    write_msp(POS_LC_insilico, "POS_LC_insilico.msp", "POS", update, profile_name)
    del POS_LC_insilico
    time.sleep(0.1)
    write_msp(POS_GC, "POS_GC.msp", "POS", update, profile_name)
    del POS_GC
    time.sleep(0.1)
    write_msp(POS_GC_insilico, "POS_GC_insilico.msp", "POS", update, profile_name)
    del POS_GC_insilico
    time.sleep(0.1)
    write_msp(NEG_LC, "NEG_LC.msp", "NEG", update, profile_name)
    del NEG_LC
    time.sleep(0.1)
    write_msp(NEG_LC_insilico, "NEG_LC_insilico.msp", "NEG", update, profile_name)
    del NEG_LC_insilico
    time.sleep(0.1)
    write_msp(NEG_GC, "NEG_GC.msp", "NEG", update, profile_name)
    del NEG_GC
    time.sleep(0.1)
    write_msp(NEG_GC_insilico, "NEG_GC_insilico.msp", "NEG", update, profile_name)
    del NEG_GC_insilico

def write_csv(df, filename, mode, update, first_run, profile_name):
    """
    :param df: pandas DataFrame object containing the data to be written to CSV
    :param filename: name of the CSV file to be created
    :param mode: mode in which the file should be opened ('w' for write, 'a' for append)
    :param update: boolean indicating whether the file should be updated or created from scratch
    :param first_run: boolean indicating whether it's the first run and headers need to be written
    :param profile_name: name of the profile associated with the CSV file
    :return: None

    This method writes a pandas DataFrame to a CSV file. The data is written in chunks, specified by the chunk_size parameter. The method uses a progress bar to track the progress of the write operation.

    The method takes the following parameters:
    - df: pandas DataFrame object containing the data to be written
    - filename: name of the CSV file to be created
    - mode: mode in which the file should be opened ('w' for write, 'a' for append)
    - update: boolean indicating whether the file should be updated or created from scratch
    - first_run: boolean indicating whether it's the first run and headers need to be written
    - profile_name: name of the profile associated with the CSV file

    The method does not return anything.
    """
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/CSV/{mode}",filename)

    chunk_size = 5000  # Taille de chaque fraction
    num_chunks = int(np.ceil(df.shape[0] / chunk_size))  # Calculer le nombre de fractions

    with tqdm(total=num_chunks, unit=" row", colour="green", desc="{:>70}".format(f"writting {filename}")) as pbar:
        for start in range(0, df.shape[0], chunk_size):
            df_slice = df[start:start + chunk_size]
            if start == 0 and first_run:
                # Write the headers for the first fraction
                df_slice.to_csv(output_file_path, mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)
            else:
                # Append to file without writing headers for other fractions
                df_slice.to_csv(output_file_path, mode='a', header=False, index=False, sep=";", quotechar='"', encoding="UTF-8")

            # Update progress bar
            pbar.update()

def writting_csv(POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df,POS_LC_df_insilico,POS_GC_df_insilico,NEG_LC_df_insilico,NEG_GC_df_insilico, first_run, profile_name, update=False):
    """
    Write dataframes to CSV files.

    :param POS_LC_df: Pandas DataFrame containing positive LC data
    :param POS_GC_df: Pandas DataFrame containing positive GC data
    :param NEG_LC_df: Pandas DataFrame containing negative LC data
    :param NEG_GC_df: Pandas DataFrame containing negative GC data
    :param POS_LC_df_insilico: Pandas DataFrame containing positive LC In Silico data
    :param POS_GC_df_insilico: Pandas DataFrame containing positive GC In Silico data
    :param NEG_LC_df_insilico: Pandas DataFrame containing negative LC In Silico data
    :param NEG_GC_df_insilico: Pandas DataFrame containing negative GC In Silico data
    :param first_run: Boolean indicating if it is the first run
    :param profile_name: Name of the profile
    :param update: Boolean indicating if the CSV files should be updated or overwritten (optional, default=False)
    :return: None
    """
    time.sleep(0.1)
    write_csv(POS_LC_df,"POS_LC.csv","POS", update, first_run, profile_name)
    del POS_LC_df
    time.sleep(0.1)
    write_csv(POS_GC_df, "POS_GC.csv","POS", update, first_run, profile_name)
    del POS_GC_df
    time.sleep(0.1)
    write_csv(NEG_LC_df, "NEG_LC.csv","NEG", update, first_run, profile_name)
    del NEG_LC_df
    time.sleep(0.1)
    write_csv(NEG_GC_df, "NEG_GC.csv","NEG", update, first_run, profile_name)
    del NEG_GC_df
    time.sleep(0.1)
    write_csv(POS_LC_df_insilico, "POS_LC_In_Silico.csv","POS", update, first_run, profile_name)
    del POS_LC_df_insilico
    time.sleep(0.1)
    write_csv(POS_GC_df_insilico, "POS_GC_In_Silico.csv","POS", update, first_run, profile_name)
    del POS_GC_df_insilico
    time.sleep(0.1)
    write_csv(NEG_LC_df_insilico, "NEG_LC_In_Silico.csv","NEG", update, first_run, profile_name)
    del NEG_LC_df_insilico
    time.sleep(0.1)
    write_csv(NEG_GC_df_insilico, "NEG_GC_In_Silico.csv","NEG", update, first_run, profile_name)
    del NEG_GC_df_insilico

def write_json(df, filename, mode, profile_name):
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
    :return: None
    """
    output_file_path = os.path.join(f"../OUTPUT/{profile_name}/JSON/{mode}", filename)

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

def writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, profile_name):
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
    time.sleep(0.1)
    write_json(POS_LC_df, "POS_LC.json", "POS", profile_name)
    del POS_LC_df
    time.sleep(0.1)
    write_json(POS_GC_df, "POS_GC.json", "POS", profile_name)
    del POS_GC_df
    time.sleep(0.1)
    write_json(NEG_LC_df, "NEG_LC.json", "NEG", profile_name)
    del NEG_LC_df
    time.sleep(0.1)
    write_json(NEG_GC_df, "NEG_GC.json", "NEG", profile_name)
    del NEG_GC_df
    time.sleep(0.1)
    write_json(POS_LC_df_insilico, "POS_LC_In_Silico.json", "POS", profile_name)
    del POS_LC_df_insilico
    time.sleep(0.1)
    write_json(POS_GC_df_insilico, "POS_GC_In_Silico.json", "POS", profile_name)
    del POS_GC_df_insilico
    time.sleep(0.1)
    write_json(NEG_LC_df_insilico, "NEG_LC_In_Silico.json", "NEG", profile_name)
    del NEG_LC_df_insilico
    time.sleep(0.1)
    write_json(NEG_GC_df_insilico, "NEG_GC_In_Silico.json", "NEG", profile_name)
    del NEG_GC_df_insilico
