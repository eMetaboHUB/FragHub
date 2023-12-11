from tqdm.auto import tqdm
from tqdm import tqdm
import pandas as pd
import time
import re

def remove_dupli_POS_LC(POS_LC):
    """
    :param POS_LC: A pandas DataFrame containing the POS_LC data.
    :return: The updated POS_LC DataFrame with duplicates removed.

    This method removes duplicate rows from the POS_LC DataFrame based on the 'INCHIKEY' and 'PEAKS_LIST' columns. It utilizes the tqdm library to display a progress bar during the removal
    * process.
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
    Remove duplicate rows from POS_LC_In_Silico DataFrame based on 'INCHIKEY' and 'PEAKS_LIST' columns.

    :param POS_LC_In_Silico: The input DataFrame containing 'INCHIKEY' and 'PEAKS_LIST' columns.
    :return: The updated DataFrame with duplicate rows removed.
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
    :param POS_GC: The input DataFrame containing the POS_GC data.
    :return: The updated DataFrame with duplicate rows removed.

    Removes duplicate rows from the input DataFrame based on the 'INCHIKEY' and 'PEAKS_LIST' columns.
    Progress bar is displayed while removing duplicates.
    Returns the updated DataFrame with duplicate rows removed.
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
    Remove duplicate rows from the POS_GC_In_Silico DataFrame.

    :param POS_GC_In_Silico: The DataFrame containing the POS_GC_In_Silico data.
    :type POS_GC_In_Silico: pandas.DataFrame
    :return: The updated DataFrame without any duplicate rows.
    :rtype: pandas.DataFrame
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
    :param NEG_LC: pandas DataFrame containing data to remove duplicates from.
    :return: pandas DataFrame with duplicates removed.

    This method removes duplicates from the given pandas DataFrame `NEG_LC` based on the columns `INCHIKEY` and `PEAKS_LIST`. It updates a progress bar using tqdm library to show the progress
    * of removing duplicates. The method returns the updated DataFrame with duplicates removed.
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
    :param NEG_LC_In_Silico: The input dataframe containing rows with duplicates to be removed.
    :return: The dataframe with duplicates removed.

    This method removes duplicate rows from the input dataframe based on two columns, 'INCHIKEY' and 'PEAKS_LIST'. It uses the tqdm library to display a progress bar during the removal process
    *. The progress bar updates as each duplicate row is removed.
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
    :param NEG_GC: pandas DataFrame containing the data to be processed.
    :return: updated NEG_GC DataFrame with duplicate rows removed.
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
    :param NEG_GC_In_Silico: A DataFrame containing data related to NEG_GC_In_Silico
    :return: The DataFrame with duplicate rows removed

    This method removes duplicate rows from the NEG_GC_In_Silico DataFrame. It takes the DataFrame as the parameter and returns the modified DataFrame with duplicates removed. The method
    * utilizes the `duplicated` method from the Pandas library to identify and remove duplicate rows based on the 'INCHIKEY' and 'PEAKS_LIST' columns. The progress of the removal process
    * is displayed using the tqdm library, which shows a progress bar indicating the number of rows processed.

    Example usage:
    ```python
    import pandas as pd
    from tqdm import tqdm

    # Create the NEG_GC_In_Silico DataFrame
    NEG_GC_In_Silico = pd.DataFrame({
        'INCHIKEY': ['key1', 'key2', 'key1', 'key3'],
        'PEAKS_LIST': ['list1', 'list2', 'list1', 'list3'],
        'VALUE': [1, 2, 3, 4]
    })

    # Remove duplicate rows
    result = remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico)
    print(result)
    ```
    ```plaintext
      INCHIKEY PEAKS_LIST  VALUE
    0     key1      list1      1
    1     key2      list2      2
    3     key3      list3      4
    ```
    """
    total_rows = len(NEG_GC_In_Silico)
    t = tqdm(total=len(NEG_GC_In_Silico), desc="NEG_GC_In_Silico", colour="green", unit=" row")

    # Supprimer les doublons et mettre à jour la barre de progression
    NEG_GC_In_Silico = NEG_GC_In_Silico.loc[~NEG_GC_In_Silico.duplicated(subset=['INCHIKEY', 'PEAKS_LIST'])]
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return NEG_GC_In_Silico

def format_comments(DF_row):
    """
    Format comments based on the given row of a DataFrame.

    :param DF_row: A row of a DataFrame.
    :type DF_row: pandas.Series
    :return: A formatted string containing information from the row.
    :rtype: str
    """
    return f'FILENAME={DF_row["FILENAME"]}; PREDICTED={DF_row["PREDICTED"]}; FRAGHUBID={DF_row["FRAGHUBID"]}; SPECTRUMID={DF_row["SPECTRUMID"]}; RESOLUTION={DF_row["RESOLUTION"]}; SYNON={DF_row["SYNON"]}; CHARGE={DF_row["CHARGE"]}; IONIZATION={DF_row["IONIZATION"]}; MSLEVEL={DF_row["MSLEVEL"]}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"]}; EXACTMASS={DF_row["EXACTMASS"]}; AVERAGEMASS={DF_row["AVERAGEMASS"]}'

def re_write_MSP_POS_LC(POS_LC_df):
    """
    :param POS_LC_df: Pandas DataFrame representing POS_LC data.
    :return: List of rewritten POS_LC data.

    This method takes a DataFrame of POS_LC data and iterates through each row. It extracts the necessary information from each row and rewrites it into a new format. The rewritten data
    * is stored in a list and returned.

    The input DataFrame, POS_LC_df, is expected to have the following columns:
    - NAME: Name of the POS_LC
    - PRECURSORMZ: Precursor m/z value
    - PRECURSORTYPE: Precursor ion type
    - FORMULA: Chemical formula
    - INCHIKEY: InChIKey
    - INCHI: InChI
    - SMILES: SMILES representation
    - RETENTIONTIME: Retention time
    - IONMODE: Ionization mode
    - INSTRUMENTTYPE: Instrument type
    - INSTRUMENT: Instrument name
    - COLLISIONENERGY: Collision energy
    - COMMENT: Additional comments
    - NUM PEAKS: Number of peaks
    - PEAKS_LIST: List of peak data

    The method uses a for loop to iterate through each row. Within the loop, it extracts the relevant information from each column and constructs a string representation of the POS_LC data
    *. The string representation is appended to the POS_LC list.

    Finally, the list of rewritten POS_LC data is returned.
    """

    POS_LC = []
    for index,row in tqdm(POS_LC_df.iterrows(), total=len(POS_LC_df), desc="\t\t  POS_LC", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC

def re_write_MSP_POS_LC_In_Silico(POS_LC_df_insilico):
    """
    Re-writes the POS_LC_in_Silico dataframe to a list of strings in the format of an MSP file.

    :param POS_LC_df_insilico: The input dataframe containing the POS_LC_in_Silico data.
    :return: A list of strings representing the re-written MSP data.
    """
    POS_LC = []
    for index, row in tqdm(POS_LC_df_insilico.iterrows(), total=len(POS_LC_df_insilico), desc="POS_LC_In_Silico", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_LC.append(SPECTRUM)

    return POS_LC


def re_write_MSP_POS_GC(POS_GC_df):
    """
    :param POS_GC_df: a pandas DataFrame containing data for POS_GC
    :return: a list of strings where each string represents a POS_GC spectrum

    The function takes a pandas DataFrame POS_GC_df as input and iterates over its rows.
    For each row, it formats the data into a string representation of a POS_GC spectrum and appends it to a list POS_GC.
    Finally, it returns the list POS_GC.
    """
    POS_GC = []
    for index, row in tqdm(POS_GC_df.iterrows(), total=len(POS_GC_df), desc="\t\t  POS_GC", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_POS_GC_In_Silico(POS_GC_df_insilico):
    """
    :param POS_GC_df_insilico: DataFrame containing in silico data for POS_GC.
    :return: List of formatted spectra.
    """
    POS_GC = []
    for index, row in tqdm(POS_GC_df_insilico.iterrows(), total=len(POS_GC_df_insilico), desc="POS_GC_In_Silico", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        POS_GC.append(SPECTRUM)

    return POS_GC

def re_write_MSP_NEG_LC(NEG_LC_df):
    """
    Re-write the data from NEG_LC_df DataFrame into a list of spectra.

    :param NEG_LC_df: The DataFrame containing the data to be re-written.
    :type NEG_LC_df: pandas.DataFrame
    :return: The list of re-written spectra.
    :rtype: list
    """
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df.iterrows(), total=len(NEG_LC_df), desc="\t\t  NEG_LC", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC

def re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico):
    """
    :param NEG_LC_df_insilico: pd.DataFrame
        The input DataFrame containing the in silico negative liquid chromatography (NEG_LC) data.
        This DataFrame should have the following columns:
        - "NAME": The name of the compound.
        - "PRECURSORMZ": The precursor m/z value.
        - "PRECURSORTYPE": The precursor type.
        - "FORMULA": The chemical formula of the compound.
        - "INCHIKEY": The InChIKey of the compound.
        - "INCHI": The InChI of the compound.
        - "SMILES": The SMILES representation of the compound.
        - "RETENTIONTIME": The retention time of the compound.
        - "IONMODE": The ion mode of the compound.
        - "INSTRUMENTTYPE": The type of instrument used.
        - "INSTRUMENT": The name of the instrument.
        - "COLLISIONENERGY": The collision energy value.
        - "NUM PEAKS": The number of peaks in the spectrum.
        - "PEAKS_LIST": A string representation of the peaks list.

    :return: list
        A list of spectra in the format specified by the given DataFrame.

    """
    NEG_LC = []
    for index, row in tqdm(NEG_LC_df_insilico.iterrows(), total=len(NEG_LC_df_insilico), desc="NEG_LC_In_Silico", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_LC.append(SPECTRUM)

    return NEG_LC


def re_write_MSP_NEG_GC(NEG_GC_df):
    """
    :param NEG_GC_df: pandas DataFrame containing the data to be processed
    :return: list of formatted spectra

    The re_write_MSP_NEG_GC method takes a pandas DataFrame, NEG_GC_df, as input. It iterates over each row of the DataFrame and formats the data into a spectrum string. The formatted spectra
    * are then added to a list, NEG_GC, which is returned as the result.

    The output spectrum string format is as follows:
    NAME: [NAME]
    PRECURSORMZ: [PRECURSORMZ]
    PRECURSORTYPE: [PRECURSORTYPE]
    FORMULA: [FORMULA]
    INCHIKEY: [INCHIKEY]
    INCHI: [INCHI]
    SMILES: [SMILES]
    RETENTIONTIME: [RETENTIONTIME]
    IONMODE: [IONMODE]
    INSTRUMENTTYPE: [INSTRUMENTTYPE]
    INSTRUMENT: [INSTRUMENT]
    COLLISIONENERGY: [COLLISIONENERGY]
    COMMENT: [COMMENTS]
    NUM PEAKS: [NUM PEAKS]
    [PEAKS_LIST]

    Example Usage:
        df = pd.DataFrame(...)
        result = re_write_MSP_NEG_GC(df)
    """
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df.iterrows(), total=len(NEG_GC_df), desc="\t\t  NEG_GC", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def re_write_MSP_NEG_GC_In_Silico(NEG_GC_df_insilico):
    """
    This method takes a DataFrame as input and converts it into a list of spectra in the MSP format.

    :param NEG_GC_df_insilico: A DataFrame containing the input data.
    :return: A list of spectra in the MSP format.
    """
    NEG_GC = []
    for index, row in tqdm(NEG_GC_df_insilico.iterrows(), total=len(NEG_GC_df_insilico), desc="NEG_GC_in_Silico", colour="green", unit=" row"):
        COMMENTS = format_comments(row)

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
        SPECTRUM = SPECTRUM + row["PEAKS_LIST"] + "\n"
        NEG_GC.append(SPECTRUM)

    return NEG_GC

def remove_duplicatas(POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico):
    """
    Remove duplicates from the provided dataframes and perform additional conversions.

    :param POS_LC: Positive LC dataframe
    :param POS_LC_In_Silico: Positive LC In Silico dataframe
    :param POS_GC: Positive GC dataframe
    :param POS_GC_In_Silico: Positive GC In Silico dataframe
    :param NEG_LC: Negative LC dataframe
    :param NEG_LC_In_Silico: Negative LC In Silico dataframe
    :param NEG_GC: Negative GC dataframe
    :param NEG_GC_In_Silico: Negative GC In Silico dataframe
    :return: A tuple containing all the cleaned dataframes.

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
