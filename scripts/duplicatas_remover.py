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

def remove_duplicatas(POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico):
    """
    Remove duplicates from the given data frames.

    :param POS_LC: a pandas DataFrame containing positive LC data
    :param POS_LC_In_Silico: a pandas DataFrame containing positive LC In Silico data
    :param POS_GC: a pandas DataFrame containing positive GC data
    :param POS_GC_In_Silico: a pandas DataFrame containing positive GC In Silico data
    :param NEG_LC: a pandas DataFrame containing negative LC data
    :param NEG_LC_In_Silico: a pandas DataFrame containing negative LC In Silico data
    :param NEG_GC: a pandas DataFrame containing negative GC data
    :param NEG_GC_In_Silico: a pandas DataFrame containing negative GC In Silico data
    :return: a tuple of pandas DataFrames with duplicates removed:
             (POS_LC_df, POS_LC_df_insilico, POS_GC_df, POS_GC_df_insilico,
             NEG_LC_df, NEG_LC_df_insilico, NEG_GC_df, NEG_GC_df_insilico)
    """
    # ========================================================================= POS_LC =========================================================================
    POS_LC_df = remove_dupli_POS_LC(POS_LC)

    # ========================================================================= POS_LC_In_Silico =========================================================================
    POS_LC_df_insilico = remove_dupli_POS_LC_In_Silico(POS_LC_In_Silico)

    # ========================================================================= POS_GC =========================================================================
    POS_GC_df = remove_dupli_POS_GC(POS_GC)

    # ========================================================================= POS_GC_In_Silico =========================================================================
    POS_GC_df_insilico = remove_dupli_POS_GC_In_Silico(POS_GC_In_Silico)

    # ========================================================================= NEG_LC =========================================================================
    NEG_LC_df = remove_dupli_NEG_LC(NEG_LC)

    # ========================================================================= NEG_LC_In_Silico =========================================================================
    NEG_LC_df_insilico = remove_dupli_NEG_LC_In_Silico(NEG_LC_In_Silico)

    # ========================================================================= NEG_GC =========================================================================
    NEG_GC_df = remove_dupli_NEG_GC(NEG_GC)

    # ========================================================================= NEG_GC_In_Silico =========================================================================
    NEG_GC_df_insilico = remove_dupli_NEG_GC_In_Silico(NEG_GC_In_Silico)


    return POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico
