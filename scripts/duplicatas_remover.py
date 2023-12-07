from tqdm.auto import tqdm
from tqdm import tqdm
import pandas as pd
import time
import re

def remove_dupli_POS_LC(POS_LC):
    """
    :param POS_LC: A pandas DataFrame containing the POS_LC data.
    :return: The POS_LC DataFrame with duplicate rows removed.

    This method takes a pandas DataFrame 'POS_LC' as input and removes duplicate rows based on the 'INCHIKEY' and 'PEAKS_LIST' columns. It updates the progress bar for each row processed
    *. The resulting DataFrame without duplicate rows is returned.
    """
    total_rows = len(POS_LC)
    t = tqdm(total=len(POS_LC), desc="\t\t  POS_LC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_LC = POS_LC.loc[~POS_LC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return  POS_LC

def remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico):
    """
    Remove duplicates from POS_LC_In_Silico DataFrame based on 'INCHIKEY' and 'PEAKS_LIST' columns.

    :param POS_LC_In_Silico: The DataFrame containing POS_LC_In_Silico data.
    :return: The DataFrame with duplicates removed.
    """
    total_rows = len(POS_LC_In_Silico)
    t = tqdm(total=len(POS_LC_In_Silico), desc="POS_LC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_LC_In_Silico = POS_LC_In_Silico.loc[~POS_LC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_LC_In_Silico

def remove_dupli_POS_GC(POS_GC):
    """
    :param POS_GC: pandas DataFrame containing the POS_GC data
    :return: a new pandas DataFrame with duplicate rows removed based on the 'INCHIKEY' and 'PEAKS_LIST' columns

    This method removes duplicate rows from the input DataFrame 'POS_GC' based on the values in the 'INCHIKEY' and 'PEAKS_LIST' columns. It updates a progress bar to track the progress of
    * the operation.

    Example usage:
        POS_GC = remove_dupli_POS_GC(POS_GC)
    """
    total_rows = len(POS_GC)
    t = tqdm(total=len(POS_GC), desc="\t\t  POS_GC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_GC = POS_GC.loc[~POS_GC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_GC

def remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico):
    """
    Removes duplicate rows from the given DataFrame `POS_GC_In_Silico` based on the 'INCHIKEY' and 'PEAKS_LIST' columns.

    :param POS_GC_In_Silico: The DataFrame containing the POS_GC_In_Silico data
    :return: The updated DataFrame with duplicate rows removed
    """
    total_rows = len(POS_GC_In_Silico)
    t = tqdm(total=len(POS_GC_In_Silico), desc="POS_GC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    POS_GC_In_Silico = POS_GC_In_Silico.loc[~POS_GC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return POS_GC_In_Silico

def remove_dupli_NEG_LC(NEG_LC):
    """
    :param NEG_LC: A pandas DataFrame containing data with columns 'INCHIKEY' and 'PEAKS_LIST'.
    :return: The pandas DataFrame with duplicate rows removed.
    """
    total_rows = len(NEG_LC)
    t = tqdm(total=len(NEG_LC), desc="\t\t  NEG_LC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_LC = NEG_LC.loc[~NEG_LC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_LC

def remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico):
    """
    Remove duplicate rows from the NEG_LC_In_Silico dataframe.

    :param NEG_LC_In_Silico: The dataframe to remove duplicates from.
    :return: The dataframe with duplicate rows removed.
    """
    total_rows = len(NEG_LC_In_Silico)
    t = tqdm(total=len(NEG_LC_In_Silico), desc="NEG_LC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_LC_In_Silico = NEG_LC_In_Silico.loc[~NEG_LC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_LC_In_Silico


def remove_dupli_NEG_GC(NEG_GC):
    """
    :param NEG_GC: DataFrame containing data with columns ['INCHIKEY', 'PEAKS_LIST']
    :return: DataFrame with duplicate rows removed based on ['INCHIKEY', 'PEAKS_LIST']

    This method removes duplicate rows from the input DataFrame based on the columns 'INCHIKEY' and 'PEAKS_LIST'.
    """
    total_rows = len(NEG_GC)
    t = tqdm(total=len(NEG_GC), desc="\t\t  NEG_GC", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_GC = NEG_GC.loc[~NEG_GC.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_GC

def remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico):
    """
    :param NEG_GC_In_Silico: A pandas DataFrame containing data with columns 'INCHIKEY' and 'PEAKS_LIST'
    :return: The modified DataFrame with duplicate rows removed

    This method removes duplicate rows from the given DataFrame based on the values in columns 'INCHIKEY' and 'PEAKS_LIST'.
    It also includes a progress bar to track the progress of the removal process.
    """
    total_rows = len(NEG_GC_In_Silico)
    t = tqdm(total=len(NEG_GC_In_Silico), desc="NEG_GC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_GC_In_Silico = NEG_GC_In_Silico.loc[~NEG_GC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_GC_In_Silico


def re_write_MSP_POS_LC(POS_LC_df):
    """
    :param POS_LC_df: DataFrame containing the POS_LC data.
    :return: A list of strings, each string representing a formatted POS_LC entry.

    This method takes in a DataFrame, `POS_LC_df`, which represents the POS_LC data. It iterates over each row in the DataFrame and formats the data into a string representation. Each string
    * contains various attributes of a POS_LC entry, separated by newline characters.

    The formatted strings are then appended to a list, `POS_LC`, and returned as the result.

    Example Usage:
        POS_LC_df = pd.read_csv('POS_LC_data.csv')  # Assuming the input data is in a CSV file
        formatted_POS_LC = re_write_MSP_POS_LC(POS_LC_df)
        for entry in formatted_POS_LC:
            print(entry)

    Output:
        FILENAME: example_filename
        PREDICTED: example_predicted
        FRAGHUBID: example_fraghubid
        ...
        NUM PEAKS: example_num_peaks
        example_peaks_list
        ...
    """
    POS_LC = []
    for index,row in tqdm(POS_LC_df.iterrows(), total=len(POS_LC_df), desc="\t\t  POS_LC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC

def re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico):
    """
    re_write_MSP_POS_LC_In_Silico method

    :param POS_LC_df_insilico: Input dataframe containing POS_LC data in silico
    :return: List of POS_LC spectrums

    This method takes an input dataframe containing POS_LC data in silico and returns a list of POS_LC spectrums. Each spectrum is generated by concatenating various attributes from the
    * input dataframe into a string representation.

    Example usage:
        POS_LC_df_insilico = pd.DataFrame(...)  # create POS_LC dataframe
        spectra = re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico)  # generate POS_LC spectrums

    """
    POS_LC = []
    for index, row in tqdm(POS_LC_df_insilico.iterrows(), total=len(POS_LC_df_insilico), desc="POS_LC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC


def re_write_MSP_POS_GC(POS_GC_df):
    """
    :param POS_GC_df: pandas DataFrame containing POS_GC data.
    :return: List of strings representing re-written POS_GC data.

    This method takes a pandas DataFrame POS_GC_df and re-writes the POS_GC data into a list of strings. Each string in the list represents a single POS_GC entry and contains various metadata
    * fields from the original data.

    Example usage:
        POS_GC_df = pd.read_csv('POS_GC.csv')
        re_written_data = re_write_MSP_POS_GC(POS_GC_df)
    """
    POS_GC = []
    for index, row in tqdm(POS_GC_df.iterrows(), total=len(POS_GC_df), desc="\t\t  POS_GC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_POS_GC_In_Silico(POS_GC_df_insilico):
    """
    Re-writes the MSP_POS_GC_In_Silico data into a list of spectra strings.

    :param POS_GC_df_insilico: DataFrame containing the MSP_POS_GC_In_Silico data
    :type POS_GC_df_insilico: pandas.DataFrame
    :return: List of spectra strings
    :rtype: list
    """
    POS_GC = []
    for index, row in tqdm(POS_GC_df_insilico.iterrows(), total=len(POS_GC_df_insilico), desc="POS_GC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_NEG_LC(NEG_LC_df):
    """
    Write a restructured text documentation for the given method `re_write_MSP_NEG_LC`.

    Parameters:
        NEG_LC_df: pandas DataFrame
            The input dataframe containing the NEG_LC data.

    Returns:
        NEG_LC: list
            A list of strings representing the rewritten MSP_NEG_LC data.

    """
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df.iterrows(), total=len(NEG_LC_df), desc="\t\t  NEG_LC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC

def re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico):
    """
    :param NEG_LC_df_insilico: Pandas DataFrame containing the in-silico negative LC data.
    :return: List of strings representing the re-written in-silico negative LC data.

    This method takes a Pandas DataFrame containing in-silico negative LC data and re-writes it in a specific format. It iterates over each row of the DataFrame and constructs a string representation
    * of the data for each row. The constructed strings are then appended to a list, which is returned at the end.

    The re-written format includes various fields from the original data, such as FILENAME, PREDICTED, FRAGHUBID, SPECTRUMID, RESOLUTION, SYNON, CHARGE, IONIZATION, MSLEVEL, FRAGMENTATION
    *MODE, NAME, PRECURSORMZ, EXACTMASS, AVERAGEMASS, PRECURSORTYPE, INSTRUMENTTYPE, INSTRUMENT, SMILES, INCHI, INCHIKEY, COLLISIONENERGY, FORMULA, RETENTIONTIME, IONMODE, COMMENT, NUM PE
    *AKS, and PEAKS_LIST. The values of these fields are concatenated to form the final string representation of each row's data.

    Note that the commented line for the PARENTMASS field is not included in the re-written format.
    """
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df_insilico.iterrows(), total=len(NEG_LC_df_insilico), desc="NEG_LC_In_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC


def re_write_MSP_NEG_GC(NEG_GC_df):
    """
    :param NEG_GC_df: dataframe containing the NEG_GC data
    :return: list of formatted NEG_GC strings
    """
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df.iterrows(), total=len(NEG_GC_df), desc="\t\t  NEG_GC", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def re_write_MSP_NEG_GC_In_Silico(NEG_GC_df_insilico):
    """
    :param NEG_GC_df_insilico: DataFrame containing in silico GC data
    :return: List of formatted spectra

    This method takes in a DataFrame containing in silico GC data and returns a list of formatted spectra. Each spectrum is formatted as a string with various attributes separated by new
    *lines. The list of spectra is stored in the NEG_GC list.
    """
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df_insilico.iterrows(), total=len(NEG_GC_df_insilico), desc="NEG_GC_in_Silico", colour="green", unit=" row"):
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGHUBID: " + row["FRAGHUBID"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        # SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
        SPECTRUM = SPECTRUM + "EXACTMASS: " + row["EXACTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "AVERAGEMASS: " + row["AVERAGEMASS"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORTYPE: " + row["PRECURSORTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENTTYPE: " + row["INSTRUMENTTYPE"] + "\n"
        SPECTRUM = SPECTRUM + "INSTRUMENT: " + row["INSTRUMENT"] + "\n"
        SPECTRUM = SPECTRUM + "SMILES: " + row["SMILES"] + "\n"
        SPECTRUM = SPECTRUM + "INCHI: " + row["INCHI"] + "\n"
        SPECTRUM = SPECTRUM + "INCHIKEY: " + row["INCHIKEY"] + "\n"
        SPECTRUM = SPECTRUM + "COLLISIONENERGY: " + row["COLLISIONENERGY"] + "\n"
        SPECTRUM = SPECTRUM + "FORMULA: " + row["FORMULA"] + "\n"
        SPECTRUM = SPECTRUM + "RETENTIONTIME: " + row["RETENTIONTIME"] + "\n"
        SPECTRUM = SPECTRUM + "IONMODE: " + row["IONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "COMMENT: " + row["COMMENT"] + "\n"
        SPECTRUM = SPECTRUM + "NUM PEAKS: " + row["NUM PEAKS"] + "\n"
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def remove_duplicatas(POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico):
    """
    :param POS_LC: the input for positive LC
    :param POS_LC_In_Silico: the input for positive LC in silico
    :param POS_GC: the input for positive GC
    :param POS_GC_In_Silico: the input for positive GC in silico
    :param NEG_LC: the input for negative LC
    :param NEG_LC_In_Silico: the input for negative LC in silico
    :param NEG_GC: the input for negative GC
    :param NEG_GC_In_Silico: the input for negative GC in silico
    :return: a tuple containing the modified POS_LC, POS_LC_df, POS_LC_df_insilico, POS_LC_In_Silico, POS_GC, POS_GC_df, POS_GC_df_insilico, POS_GC_In_Silico, NEG_LC, NEG_LC_df, NEG_LC_df
    *_insilico, NEG_LC_In_Silico, NEG_GC, NEG_GC_df, NEG_GC_df_insilico, NEG_GC_In_Silico


    """
    # ========================================================================= POS_LC =========================================================================
    POS_LC_df = remove_dupli_POS_LC(POS_LC)
    # Re convert to MSP
    POS_LC = re_write_MSP_POS_LC(POS_LC_df)

    # ========================================================================= POS_LC_In_Silico =========================================================================
    POS_LC_df_insilico = remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico)
    # Re convert to MSP
    POS_LC_In_Silico = re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico)

    # ========================================================================= POS_GC =========================================================================
    POS_GC_df = remove_dupli_POS_GC(POS_GC)
    # Re convert to MSP
    POS_GC = re_write_MSP_POS_GC(POS_GC_df)

    # ========================================================================= POS_GC_In_Silico =========================================================================
    POS_GC_df_insilico = remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico)
    # Re convert to MSP
    POS_GC_In_Silico = re_write_MSP_POS_GC_In_Silico(POS_GC_df_insilico)

    # ========================================================================= NEG_LC =========================================================================
    NEG_LC_df = remove_dupli_NEG_LC(NEG_LC)
    # Re convert to MSP
    NEG_LC = re_write_MSP_NEG_LC(NEG_LC_df)

    # ========================================================================= NEG_LC_In_Silico =========================================================================
    NEG_LC_df_insilico = remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico)
    # Re convert to MSP
    NEG_LC_In_Silico = re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico)

    # ========================================================================= NEG_GC =========================================================================
    NEG_GC_df = remove_dupli_NEG_GC(NEG_GC)
    # Re convert to MSP
    NEG_GC = re_write_MSP_NEG_GC(NEG_GC_df)

    # ========================================================================= NEG_GC_In_Silico =========================================================================
    NEG_GC_df_insilico = remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico)
    # Re convert to MSP
    NEG_GC_In_Silico = re_write_MSP_NEG_GC_In_Silico(NEG_GC_df_insilico)


    return POS_LC,POS_LC_df,POS_LC_df_insilico,POS_LC_In_Silico,POS_GC,POS_GC_df,POS_GC_df_insilico,POS_GC_In_Silico,NEG_LC,NEG_LC_df,NEG_LC_df_insilico,NEG_LC_In_Silico,NEG_GC,NEG_GC_df,NEG_GC_df_insilico,NEG_GC_In_Silico
