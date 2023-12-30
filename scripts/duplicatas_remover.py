from tqdm.auto import tqdm
from tqdm import tqdm
import pandas as pd
import time
import re

def removing_duplicates(dataframe, name):
    """
    Remove duplicates from a dataframe and update the progress bar.

    :param dataframe: The dataframe to remove duplicates from.
    :param name: The name to display on the progress bar.
    :return: The dataframe with duplicates removed.
    """
    total_rows = len(dataframe)
    t = tqdm(total=len(dataframe), desc=name, colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    dataframe = dataframe.loc[~dataframe.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return dataframe

def remove_duplicatas(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df):
    """
    Remove duplicates from the given dataframes.

    :param POS_LC_df: DataFrame containing positive LC data.
    :param POS_LC_In_Silico_df: DataFrame containing positive LC In Silico data.
    :param POS_GC_df: DataFrame containing positive GC data.
    :param POS_GC_In_Silico_df: DataFrame containing positive GC In Silico data.
    :param NEG_LC_df: DataFrame containing negative LC data.
    :param NEG_LC_In_Silico_df: DataFrame containing negative LC In Silico data.
    :param NEG_GC_df: DataFrame containing negative GC data.
    :param NEG_GC_In_Silico_df: DataFrame containing negative GC In Silico data.
    :return: Tuple containing the updated dataframes without duplicates.
    """
    # ========================================================================= POS_LC =========================================================================
    POS_LC_df = removing_duplicates(POS_LC_df, "POS_LC")

    # ========================================================================= POS_LC_In_Silico =========================================================================
    POS_LC_In_Silico_df = removing_duplicates(POS_LC_In_Silico_df, "POS_LC_In_Silico")

    # ========================================================================= POS_GC =========================================================================
    POS_GC_df = removing_duplicates(POS_GC_df, "POS_GC")

    # ========================================================================= POS_GC_In_Silico =========================================================================
    POS_GC_In_Silico_df = removing_duplicates(POS_GC_In_Silico_df, "POS_GC_In_Silico")

    # ========================================================================= NEG_LC =========================================================================
    NEG_LC_df = removing_duplicates(NEG_LC_df, "NEG_LC")

    # ========================================================================= NEG_LC_In_Silico =========================================================================
    NEG_LC_In_Silico_df = removing_duplicates(NEG_LC_In_Silico_df, "NEG_LC_In_Silico")

    # ========================================================================= NEG_GC =========================================================================
    NEG_GC_df = removing_duplicates(NEG_GC_df, "NEG_GC")

    # ========================================================================= NEG_GC_In_Silico =========================================================================
    NEG_GC_In_Silico_df = removing_duplicates(NEG_GC_In_Silico_df, "NEG_GC_In_Silico")


    return POS_LC_df,POS_LC_In_Silico_df,POS_GC_df,POS_GC_In_Silico_df,NEG_LC_df,NEG_LC_In_Silico_df,NEG_GC_df,NEG_GC_In_Silico_df
