from tqdm.auto import tqdm
from tqdm import tqdm
import pandas as pd
import os

def removing_duplicates(dataframe, name, mode, update, first_run, profile_name):
    """
    Remove duplicates from a dataframe and update the progress bar.

    :param dataframe: The input dataframe.
    :param name: The name of the dataframe.
    :param mode: The mode of the dataframe.
    :param update: A boolean indicating whether to update the dataframe.
    :param first_run: A boolean indicating whether it is the first run.
    :param profile_name: The profile name.
    :return: The dataframe after removing duplicates.
    """
    if update and not first_run:
        filename = os.path.abspath(f"../OUTPUT/{profile_name}/CSV/{mode}/{name}.csv")
        if os.path.exists(filename):
            csv_dataframe = pd.read_csv(filename, sep=";", quotechar='"', encoding="UTF-8")
        else:
            csv_dataframe = pd.DataFrame()

        dataframe = pd.concat([dataframe, csv_dataframe])

    total_rows = len(dataframe)
    t = tqdm(total=len(dataframe), desc="{:>70}".format(name), colour="green", unit=" row")

    # Delete duplicates and update progress bar
    dataframe = dataframe.loc[~dataframe.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Close progress bar
    t.close()

    return dataframe

def remove_duplicatas(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df, first_run, profile_name, update=False):
    """
    Remove duplicates from the given dataframes.

    :param POS_LC_df: The dataframe for positive LC data.
    :param POS_LC_In_Silico_df: The dataframe for positive LC In Silico data.
    :param POS_GC_df: The dataframe for positive GC data.
    :param POS_GC_In_Silico_df: The dataframe for positive GC In Silico data.
    :param NEG_LC_df: The dataframe for negative LC data.
    :param NEG_LC_In_Silico_df: The dataframe for negative LC In Silico data.
    :param NEG_GC_df: The dataframe for negative GC data.
    :param NEG_GC_In_Silico_df: The dataframe for negative GC In Silico data.
    :param first_run: A boolean value indicating if it's the first run.
    :param profile_name: The name of the profile.
    :param update: A boolean value indicating if the dataframes should be updated. Default is False.
    :return: The updated dataframes for POS_LC, POS_LC_In_Silico, POS_GC, POS_GC_In_Silico, NEG_LC, NEG_LC_In_Silico, NEG_GC, and NEG_GC_In_Silico.
    """
    # ========================================================================= POS_LC =========================================================================
    POS_LC_df = removing_duplicates(POS_LC_df, "POS_LC", "POS", update, first_run, profile_name)

    # ========================================================================= POS_LC_In_Silico =========================================================================
    POS_LC_In_Silico_df = removing_duplicates(POS_LC_In_Silico_df, "POS_LC_In_Silico", "POS", update, first_run, profile_name)

    # ========================================================================= POS_GC =========================================================================
    POS_GC_df = removing_duplicates(POS_GC_df, "POS_GC", "POS", update, first_run, profile_name)

    # ========================================================================= POS_GC_In_Silico =========================================================================
    POS_GC_In_Silico_df = removing_duplicates(POS_GC_In_Silico_df, "POS_GC_In_Silico", "POS", update, first_run, profile_name)

    # ========================================================================= NEG_LC =========================================================================
    NEG_LC_df = removing_duplicates(NEG_LC_df, "NEG_LC", "NEG", update, first_run, profile_name)

    # ========================================================================= NEG_LC_In_Silico =========================================================================
    NEG_LC_In_Silico_df = removing_duplicates(NEG_LC_In_Silico_df, "NEG_LC_In_Silico", "NEG", update, first_run, profile_name)

    # ========================================================================= NEG_GC =========================================================================
    NEG_GC_df = removing_duplicates(NEG_GC_df, "NEG_GC", "NEG", update, first_run, profile_name)

    # ========================================================================= NEG_GC_In_Silico =========================================================================
    NEG_GC_In_Silico_df = removing_duplicates(NEG_GC_In_Silico_df, "NEG_GC_In_Silico", "NEG", update, first_run, profile_name)


    return POS_LC_df,POS_LC_In_Silico_df,POS_GC_df,POS_GC_In_Silico_df,NEG_LC_df,NEG_LC_In_Silico_df,NEG_GC_df,NEG_GC_In_Silico_df
