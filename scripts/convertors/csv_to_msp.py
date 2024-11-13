from tqdm import tqdm
import pandas as pd
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
    # If it exists, include the attribute value. If not, include 'UNKNOWN'.
    return f'FILENAME={DF_row["FILENAME"] if DF_row["FILENAME"] else "UNKNOWN"}; PREDICTED={DF_row["PREDICTED"] if DF_row["PREDICTED"] else "UNKNOWN"}; SPLASH={DF_row["SPLASH"] if DF_row["SPLASH"] else "UNKNOWN"}; SPECTRUMID={DF_row["SPECTRUMID"] if DF_row["SPECTRUMID"] else "UNKNOWN"}; RESOLUTION={DF_row["RESOLUTION"] if DF_row["RESOLUTION"] else "UNKNOWN"}; SYNON={DF_row["SYNON"] if DF_row["SYNON"] else "UNKNOWN"}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"] if DF_row["FRAGMENTATIONMODE"] else "UNKNOWN"}; AVERAGEMASS={DF_row["AVERAGEMASS"] if DF_row["AVERAGEMASS"] else "UNKNOWN"}; ENTROPY={DF_row["ENTROPY"] if DF_row["ENTROPY"] else "UNKNOWN"}'

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
    # Convert the data type of all columns in the dataframe to string
    dataframe = dataframe.astype(str)

    # Initialize an empty list to store the formatted spectra strings
    spectrum_list = []

    # Iterate over each row in the dataframe
    with tqdm(total=len(dataframe), desc="{:>70}".format(name), colour="green", unit=" row") as pbar:
        for index, row in dataframe.iterrows():
            # Format the comments for the current row
            COMMENTS = format_comments(row)

            # Here begins the creation of the spectrum string based on the values in the row
            # Each attribute is checked if it exists in the row; if it doesn't, the value "UNKNOWN" is used
            # The following block does this for several attributes such as NAME, PRECURSORMZ, PRECURSORTYPE, etc.
            SPECTRUM = ""
            SPECTRUM = SPECTRUM + "NAME: " + (row["NAME"] if row["NAME"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORMZ: " + (row["PRECURSORMZ"] if row["PRECURSORMZ"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + (row["PRECURSORTYPE"] if row["PRECURSORTYPE"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "FORMULA: " + (row["FORMULA"] if row["FORMULA"] else "UNKNOWN") + "\n"
            # Ontology ???
            SPECTRUM = SPECTRUM + "INCHIKEY: " + (row["INCHIKEY"] if row["INCHIKEY"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "INCHI: " + (row["INCHI"] if row["INCHI"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "SMILES: " + (row["SMILES"] if row["SMILES"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "RETENTIONTIME: " + (row["RETENTIONTIME"] if row["RETENTIONTIME"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "IONMODE: " + (row["IONMODE"] if row["IONMODE"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + (row["INSTRUMENTTYPE"] if row["INSTRUMENTTYPE"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "INSTRUMENT: " + (row["INSTRUMENT"] if row["INSTRUMENT"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + (row["COLLISIONENERGY"] if row["COLLISIONENERGY"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "EXACTMASS: " + (str(row["EXACTMASS"]) if row["EXACTMASS"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "IONIZATION: " + (row["IONIZATION"] if row["IONIZATION"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "MSLEVEL: " + (str(row["MSLEVEL"]) if row["MSLEVEL"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "COMMENT: " + (COMMENTS if COMMENTS else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + "NUM PEAKS: " + (row["NUM PEAKS"] if row["NUM PEAKS"] else "UNKNOWN") + "\n"
            SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
            # Finally, append the formatted spectrum string for the current row to the spectrum_list
            spectrum_list.append(SPECTRUM)
            pbar.update()
        pbar.close()

    # Return the list of formatted spectra strings
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
    # Create MSP file from POS_LC dataframe
    POS_LC = dataframe_to_msp(POS_LC_df, "POS_LC")
    # Create MSP file from POS_LC insilico dataframe
    POS_LC_insilico = dataframe_to_msp(POS_LC_df_insilico, "POS_LC_insilico")
    # Create MSP file from POS_GC dataframe
    POS_GC = dataframe_to_msp(POS_GC_df, "POS_GC")
    # Create MSP file from POS_GC insilico dataframe
    POS_GC_insilico = dataframe_to_msp(POS_GC_df_insilico, "POS_GC_insilico")
    # Create MSP file from NEG_LC dataframe
    NEG_LC = dataframe_to_msp(NEG_LC_df, "NEG_LC")
    # Create MSP file from NEG_LC insilico dataframe
    NEG_LC_insilico = dataframe_to_msp(NEG_LC_df_insilico, "NEG_LC_insilico")
    # Create MSP file from NEG_GC dataframe
    NEG_GC = dataframe_to_msp(NEG_GC_df, "NEG_GC")
    # Create MSP file from NEG_GC insilico dataframe
    NEG_GC_insilico = dataframe_to_msp(NEG_GC_df_insilico, "NEG_GC_insilico")

    # Return original dataframes as well as the generated MSP files.
    return POS_LC_df, POS_LC, POS_LC_df_insilico, POS_LC_insilico, POS_GC_df, POS_GC, POS_GC_df_insilico, POS_GC_insilico, NEG_LC_df, NEG_LC, NEG_LC_df_insilico, NEG_LC_insilico, NEG_GC_df, NEG_GC, NEG_GC_df_insilico, NEG_GC_insilico
