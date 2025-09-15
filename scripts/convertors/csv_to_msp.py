import pandas as pd
import time
import re

def format_comments(DF_row):
    """
    Format the comments from a DataFrame row.

    :param DF_row: A row from a DataFrame.
    :type DF_row: pandas.Series
    :return: The formatted comments string.
    :rtype: str
    """

    # Use a formatted string (f-string) to create the output string.
    # For each attribute, use a ternary conditional operator to check if the attribute value exists in the row.
    # If it exists, include the attribute value. If not, include 'NOT FOUND'.
    return f'FILENAME={DF_row["FILENAME"] if DF_row["FILENAME"] else "NOT FOUND"}; FILEHASH={DF_row["FILEHASH"] if DF_row["FILEHASH"] else "NOT FOUND"}; PREDICTED={DF_row["PREDICTED"] if DF_row["PREDICTED"] else "NOT FOUND"}; SPLASH={DF_row["SPLASH"] if DF_row["SPLASH"] else "NOT FOUND"}; SPECTRUMID={DF_row["SPECTRUMID"] if DF_row["SPECTRUMID"] else "NOT FOUND"}; RESOLUTION={DF_row["RESOLUTION"] if DF_row["RESOLUTION"] else "NOT FOUND"}; SYNON={DF_row["SYNON"] if DF_row["SYNON"] else "NOT FOUND"}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"] if DF_row["FRAGMENTATIONMODE"] else "NOT FOUND"}; AVERAGEMASS={DF_row["AVERAGEMASS"] if DF_row["AVERAGEMASS"] else "NOT FOUND"}; ENTROPY={DF_row["ENTROPY"] if DF_row["ENTROPY"] else "NOT FOUND"}; ONTOLOGIES = "CLASSYFIRE_SUPERCLASS={DF_row['CLASSYFIRE_SUPERCLASS'] if DF_row['CLASSYFIRE_SUPERCLASS'] else 'NOT FOUND'}, CLASSYFIRE_CLASS = {DF_row['CLASSYFIRE_CLASS'] if DF_row['CLASSYFIRE_CLASS'] else 'NOT FOUND'}, CLASSYFIRE_SUBCLASS = {DF_row['CLASSYFIRE_SUBCLASS'] if DF_row['CLASSYFIRE_SUBCLASS'] else 'NOT FOUND'}, NPCLASS_PATHWAY = {DF_row['NPCLASS_PATHWAY'] if DF_row['NPCLASS_PATHWAY'] else 'NOT FOUND'}, NPCLASS_SUPERCLASS = {DF_row['NPCLASS_SUPERCLASS'] if DF_row['NPCLASS_SUPERCLASS'] else 'NOT FOUND'}, NPCLASS_CLASS = {DF_row['NPCLASS_CLASS'] if DF_row['NPCLASS_CLASS'] else 'NOT FOUND'}"'

def dataframe_to_msp(dataframe, name, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Converts a DataFrame containing spectral data into a list of formatted strings representing the spectra.

    :param dataframe: A pandas DataFrame containing spectral data.
    :param name: The name of the spectrum (used as context for the callbacks).
    :param progress_callback: A callable function to notify progress.
    :param total_items_callback: A callable function to set the total number of items.
    :param prefix_callback: A callable function to signal the context of the task.
    :param item_type_callback: A callable function to signal the type of items being processed.
    :return: A list of formatted strings representing the spectra.
    """
    # Reindex the DataFrame to ensure consecutive indexing (important after filtering/splitting)
    dataframe = dataframe.reset_index(drop=True)

    # Convert all columns to string type
    dataframe = dataframe.astype(str)

    # Indicate the ongoing task
    if prefix_callback:
        prefix_callback(f"Formatting {name} to MSP:")

    # Define the type of items being processed
    if item_type_callback:
        item_type_callback("rows")

    # Signal the total number of items to process
    total_rows = len(dataframe)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Initialize completed items

    # Initialize a list for the formatted strings
    spectrum_list = []

    # Iterate over each row of the DataFrame
    for index, row in dataframe.iterrows():
        # Format the comment for the current row
        COMMENTS = format_comments(row)

        # Create the spectrum string based on the row values
        SPECTRUM = ""
        SPECTRUM += "NAME: " + (row["NAME"] if row["NAME"] else "NOT FOUND") + "\n"
        SPECTRUM += "PRECURSORMZ: " + (row["PRECURSORMZ"] if row["PRECURSORMZ"] else "NOT FOUND") + "\n"
        SPECTRUM += "PRECURSORTYPE: " + (row["PRECURSORTYPE"] if row["PRECURSORTYPE"] else "NOT FOUND") + "\n"
        SPECTRUM += "FORMULA: " + (row["FORMULA"] if row["FORMULA"] else "NOT FOUND") + "\n"
        SPECTRUM += "INCHIKEY: " + (row["INCHIKEY"] if row["INCHIKEY"] else "NOT FOUND") + "\n"
        SPECTRUM += "INCHI: " + (row["INCHI"] if row["INCHI"] else "NOT FOUND") + "\n"
        SPECTRUM += "SMILES: " + (row["SMILES"] if row["SMILES"] else "NOT FOUND") + "\n"
        SPECTRUM += "RT: " + (row["RT"] if row["RT"] else "NOT FOUND") + "\n"
        SPECTRUM += "IONMODE: " + (row["IONMODE"] if row["IONMODE"] else "NOT FOUND") + "\n"
        SPECTRUM += "INSTRUMENTTYPE: " + (row["INSTRUMENTTYPE"] if row["INSTRUMENTTYPE"] else "NOT FOUND") + "\n"
        SPECTRUM += "INSTRUMENT: " + (row["INSTRUMENT"] if row["INSTRUMENT"] else "NOT FOUND") + "\n"
        SPECTRUM += "COLLISIONENERGY: " + (row["COLLISIONENERGY"] if row["COLLISIONENERGY"] else "NOT FOUND") + "\n"
        SPECTRUM += "EXACTMASS: " + (row["EXACTMASS"] if row["EXACTMASS"] else "NOT FOUND") + "\n"
        SPECTRUM += "IONIZATION: " + (row["IONIZATION"] if row["IONIZATION"] else "NOT FOUND") + "\n"
        SPECTRUM += "MSLEVEL: " + (row["MSLEVEL"] if row["MSLEVEL"] else "NOT FOUND") + "\n"
        SPECTRUM += "COMMENT: " + (COMMENTS if COMMENTS else "NOT FOUND") + "\n"
        SPECTRUM += "NUM PEAKS: " + (row["NUM PEAKS"] if row["NUM PEAKS"] else "NOT FOUND") + "\n"
        SPECTRUM += row["PEAKS_LIST"] + "\n"

        # Add the formatted spectrum to the list of spectra
        spectrum_list.append(SPECTRUM)

        # Update the progress
        if progress_callback:
            progress_callback(index + 1)

    # Return the list of formatted spectra
    return spectrum_list

def csv_to_msp(POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    :param POS_LC_df: DataFrame containing positive LC data
    :param POS_LC_df_insilico: DataFrame containing positive LC insilico data
    :param POS_GC_df: DataFrame containing positive GC data
    :param POS_GC_df_insilico: DataFrame containing positive GC insilico data
    :param NEG_LC_df: DataFrame containing negative LC data
    :param NEG_LC_df_insilico: DataFrame containing negative LC insilico data
    :param NEG_GC_df: DataFrame containing negative GC data
    :param NEG_GC_df_insilico: DataFrame containing negative GC insilico data
    :return: Tuple containing the following dataframes and MSP files:
        - POS_LC_df: Positive LC DataFrame
        - POS_LC: Positive LC MSP file
        - POS_LC_df_insilico: Positive LC insilico DataFrame
        - POS_LC_insilico: Positive LC insilico MSP file
        - POS_GC_df: Positive GC DataFrame
        - POS_GC: Positive GC MSP file
        - POS_GC_df_insilico: Positive GC insilico DataFrame
        - POS_GC_insilico: Positive GC insilico MSP file
        - NEG_LC_df: Negative LC DataFrame
        - NEG_LC: Negative LC MSP file
        - NEG_LC_df_insilico: Negative LC insilico DataFrame
        - NEG_LC_insilico: Negative LC insilico MSP file
        - NEG_GC_df: Negative GC DataFrame
        - NEG_GC: Negative GC MSP file
        - NEG_GC_df_insilico: Negative GC insilico DataFrame
        - NEG_GC_insilico: Negative GC insilico MSP file

    """
    # Create MSP file from POS_LC dataframe
    time.sleep(0.1)
    POS_LC = dataframe_to_msp(POS_LC_df, "POS_LC", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from POS_LC insilico dataframe
    time.sleep(0.1)
    POS_LC_insilico = dataframe_to_msp(POS_LC_df_insilico, "POS_LC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from POS_GC dataframe
    time.sleep(0.1)
    POS_GC = dataframe_to_msp(POS_GC_df, "POS_GC", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from POS_GC insilico dataframe
    time.sleep(0.1)
    POS_GC_insilico = dataframe_to_msp(POS_GC_df_insilico, "POS_GC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from NEG_LC dataframe
    time.sleep(0.1)
    NEG_LC = dataframe_to_msp(NEG_LC_df, "NEG_LC", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from NEG_LC insilico dataframe
    time.sleep(0.1)
    NEG_LC_insilico = dataframe_to_msp(NEG_LC_df_insilico, "NEG_LC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from NEG_GC dataframe
    time.sleep(0.1)
    NEG_GC = dataframe_to_msp(NEG_GC_df, "NEG_GC", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Create MSP file from NEG_GC insilico dataframe
    time.sleep(0.1)
    NEG_GC_insilico = dataframe_to_msp(NEG_GC_df_insilico, "NEG_GC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Return original dataframes as well as the generated MSP files.
    return POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico
