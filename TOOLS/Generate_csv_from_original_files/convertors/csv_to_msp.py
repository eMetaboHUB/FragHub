from tqdm import tqdm
import pandas as pd
import ast
import re

def format_comments(DF_row):
    """
    Format comments based on the given row of a DataFrame.

    :param DF_row: A row of a DataFrame.
    :type DF_row: pandas.Series
    :return: A formatted string containing information from the row.
    :rtype: str
    """
    return f'FILENAME={DF_row["FILENAME"]}; PREDICTED={DF_row["PREDICTED"]}; FRAGHUBID={DF_row["FRAGHUBID"]}; SPECTRUMID={DF_row["SPECTRUMID"]}; RESOLUTION={DF_row["RESOLUTION"]}; SYNON={DF_row["SYNON"]}; IONIZATION={DF_row["IONIZATION"]}; MSLEVEL={DF_row["MSLEVEL"]}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"]}; EXACTMASS={DF_row["EXACTMASS"]}; AVERAGEMASS={DF_row["AVERAGEMASS"]}'

def dataframe_to_msp(dataframe, name):
    """
    :param dataframe: A pandas DataFrame containing spectral data.
    :param name: The name of the spectrum. Used for tqdm progress bar description.
    :return: A list of formatted spectra strings.

    This method takes a DataFrame containing spectral data and converts it into a list of formatted spectra strings. Each row in the DataFrame represents a spectrum, and the columns contain
    * different attributes of the spectrum such as name, precursor m/z, precursor type, etc.

    The method iterates over each row in the DataFrame and formats the spectrum information into a string using a predefined format. The resulting formatted spectrum strings are appended
    * to the spectrum_list.

    The method uses the tqdm library to display a progress bar while iterating over the rows of the DataFrame. The name parameter is used as the description for the progress bar.

    Example usage:
        dataframe = pd.DataFrame(...)  # Create a DataFrame with spectral data
        name = "Spectrum Conversion"  # Specify a name for the spectrum
        spectra = dataframe_to_msp(dataframe, name)  # Convert DataFrame to list of formatted spectra strings
    """
    spectrum_list = []
    for index, row in tqdm(dataframe.iterrows(), total=len(dataframe), desc="{:>70}".format(name), colour="green", unit=" row"):
        COMMENTS = format_comments(row)
        peak_list = ast.literal_eval(row["PEAKS_LIST"])

        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        # Ontology ???
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + COMMENTS + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + "\n".join(["\t".join(map(str, sub_list)) for sub_list in peak_list]) + "\n"
        spectrum_list.append(SPECTRUM)

    return spectrum_list

def csv_to_msp(POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico):
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
    POS_LC = dataframe_to_msp(POS_LC_df, "POS_LC")
    POS_LC_insilico = dataframe_to_msp(POS_LC_df_insilico, "POS_LC_insilico")
    POS_GC = dataframe_to_msp(POS_GC_df, "POS_GC")
    POS_GC_insilico = dataframe_to_msp(POS_GC_df_insilico, "POS_GC_insilico")
    NEG_LC = dataframe_to_msp(NEG_LC_df, "NEG_LC")
    NEG_LC_insilico = dataframe_to_msp(NEG_LC_df_insilico, "NEG_LC_insilico")
    NEG_GC = dataframe_to_msp(NEG_GC_df, "NEG_GC")
    NEG_GC_insilico = dataframe_to_msp(NEG_GC_df_insilico, "NEG_GC_insilico")


    return POS_LC_df,POS_LC,POS_LC_df_insilico,POS_LC_insilico,POS_GC_df,POS_GC,POS_GC_df_insilico,POS_GC_insilico,NEG_LC_df,NEG_LC,NEG_LC_df_insilico,NEG_LC_insilico,NEG_GC_df,NEG_GC,NEG_GC_df_insilico,NEG_GC_insilico
