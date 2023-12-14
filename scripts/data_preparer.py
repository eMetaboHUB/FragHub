import concurrent.futures
from tqdm import tqdm
import pandas as pd
import time


def peak_list_to_float(spectrum_df):
    spectrum_df.at[0, "peak_list"] = spectrum_df.at[0, "peak_list"].astype(float)
    return spectrum_df

def prepare_data(spectrum_df):
    time.sleep(0.000000001)  # Needed to ensure progress bar display update (1ns)
    spectrum_df = peak_list_to_float(spectrum_df)

    return spectrum_df

def data_preparer_processing(spectrum_df_list):
    """
    :param spectrum_df_list: A list of spectrum dataframes to be processed.
    :return: A list containing the results of the processing.
    """


    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(prepare_data, spectrum_df_list), total=len(spectrum_df_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.