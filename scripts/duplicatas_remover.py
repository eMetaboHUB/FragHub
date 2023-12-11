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
    return f'FILENAME={DF_row["FILENAME"]}; PREDICTED={DF_row["PREDICTED"]}; FRAGHUBID={DF_row["FRAGHUBID"]}; SPECTRUMID={DF_row["SPECTRUMID"]}; RESOLUTION={DF_row["RESOLUTION"]}; SYNON={DF_row["SYNON"]}; CHARGE={DF_row["CHARGE"]}; IONIZATION={DF_row["IONIZATION"]}; MSLEVEL={DF_row["MSLEVEL"]}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"]}; EXACTMASS={DF_row["EXACTMASS"]}; AVERAGEMASS={DF_row["AVERAGEMASS"]}'

def re_write_MSP_POS_LC(POS_LC_df):
    """
    Re-writes the MSP_POS_LC DataFrame to a list of strings.

    :param POS_LC_df: DataFrame containing the MSP_POS_LC data.
    :return: A list of strings representing the re-written POS_LC data.
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
    re_write_MSP_POS_LC_In_Silico

    :param POS_LC_df_insilico: pandas DataFrame containing the POS_LC insilico data
    :return: list of modified POS_LC insilico data

    This method takes in a pandas DataFrame containing POS_LC insilico data and modifies it by concatenating specific columns into a single string, then appends the resulting string to a
    * list. The resulting list is then returned.

    Example usage:
    ---------------
    import pandas as pd

    # Create a sample DataFrame
    data = {
        "FILENAME": ["file1", "file2"],
        "PREDICTED": ["predicted1", "predicted2"],
        "FRAGHUBID": ["id1", "id2"],
        # Add other columns here
    }

    df = pd.DataFrame(data)

    # Call the method and store the result
    modified_data = re_write_MSP_POS_LC_In_Silico(df)
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
    :param POS_GC_df: DataFrame containing POS_GC data
    :return: List of POS_GC strings

    This method takes a DataFrame `POS_GC_df` as input and iterates over each row in the DataFrame. For each row, it constructs a POS_GC string and appends it to a list. The constructed
    * POS_GC strings include various information from the row's columns, separated by newline characters. The resulting list of POS_GC strings is returned as the output of the method.
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
    :param POS_GC_df_insilico: DataFrame containing the POS_GC data.
    :return: List of strings containing the transformed POS_GC data.

    This method takes a DataFrame POS_GC_df_insilico as input and transforms it into a list of strings representing the POS_GC data. Each row in the DataFrame is processed to create a string
    * representation in the SPECTRUM format. The resulting strings are appended to a list called POS_GC, which is then returned.
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
    :param NEG_LC_df: pandas DataFrame containing the data for NEG_LC
    :return: List of strings representing the re-written SPECTRUM data

    The `re_write_MSP_NEG_LC` function takes a pandas DataFrame `NEG_LC_df` as input and returns a list of strings `NEG_LC` representing the re-written SPECTRUM data.

    The function iterates over each row in `NEG_LC_df` using the `iterrows()` method. For each row, it constructs a string `SPECTRUM` by concatenating various columns of the row with relevant
    * labels. The constructed `SPECTRUM` string is then appended to the `NEG_LC` list.

    Finally, the function returns the `NEG_LC` list containing all the re-written SPECTRUM data.
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
    :param NEG_LC_df_insilico: A pandas DataFrame containing in silico negative LC data.
    :return: A list of formatted spectra.

    The re_write_MSP_NEG_LC_In_Silico method takes a DataFrame of in silico negative LC data and returns a list of formatted spectra.

    Each row in the input DataFrame represents a single spectrum and contains various columns representing different attributes of the spectrum. The method iterates over each row of the
    * DataFrame and formats the spectrum information into a string, following a specific pattern. The formatted string is then added to a list.

    The formatted string for each spectrum includes the following information:
    - FILENAME
    - PREDICTED
    - FRAGHUBID
    - SPECTRUMID
    - RESOLUTION
    - SYNON
    - CHARGE
    - IONIZATION
    - MSLEVEL
    - FRAGMENTATIONMODE
    - NAME
    - PRECURSORMZ
    - EXACTMASS
    - AVERAGEMASS
    - PRECURSORTYPE
    - INSTRUMENTTYPE
    - INSTRUMENT
    - SMILES
    - INCHI
    - INCHIKEY
    - COLLISIONENERGY
    - FORMULA
    - RETENTIONTIME
    - IONMODE
    - COMMENT
    - NUM PEAKS
    - PEAKS_LIST

    The list of formatted spectra is then returned.

    Example usage:
    ```
    import pandas as pd

    # Create a DataFrame containing the in silico negative LC data
    NEG_LC_df_insilico = pd.DataFrame({
        "FILENAME": ["spectrum1", "spectrum2"],
        "PREDICTED": ["predicted1", "predicted2"],
        "FRAGHUBID": ["fraghubid1", "fraghubid2"],
        ...
        "PEAKS_LIST": ["peaks_list1", "peaks_list2"]
    })

    # Call the re_write_MSP_NEG_LC_In_Silico method
    formatted_spectra = re_write_MSP_NEG_LC_In_Silico(NEG_LC_df_insilico)
    ```
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
    :param NEG_GC_df: pandas dataframe containing the negative global compounds (NEG_GC) data
    :return: list of strings representing the rewritten NEG_GC data

    This method takes a pandas dataframe containing the negative global compounds (NEG_GC) data and rewrites it in a specific format. It iterates through each row of the dataframe, concaten
    *ates the row values into a string called SPECTRUM, and appends it to the NEG_GC list. Finally, it returns the updated NEG_GC list.

    Note: The tqdm function is used to provide a progress bar for the iteration process.
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
    :param NEG_GC_df_insilico: DataFrame containing in silico negative GC spectra data
    :return: List of formatted spectra strings

    This method takes a DataFrame containing in silico negative GC spectra data and generates a formatted spectrum string for each row in the DataFrame. The formatted spectrum strings are
    * appended to a list and returned.

    The method iterates over each row in the DataFrame using the `iterrows()` method. For each row, it creates a new spectrum string `SPECTRUM` and appends various properties from the row
    * to the string using concatenation. Once all properties are appended, the spectrum string is added to the `NEG_GC` list.

    Finally, the method returns the `NEG_GC` list containing all the formatted spectra strings.
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
