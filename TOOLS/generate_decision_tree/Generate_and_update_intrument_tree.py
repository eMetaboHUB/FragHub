import pandas as pd
import itertools
import copy

def gen_missing_comb_multi(df):
    """
    Generate missing combinations for multiple columns.

    :param df: The DataFrame containing the columns.
    :return: A new DataFrame with additional rows representing the missing combinations.

    """
    new_rows = []
    columns = df.columns.tolist()
    for col_num in range(1, len(columns) + 1):
        for col_comb in itertools.combinations(columns, col_num):
            for i, row in df.iterrows():
                new_row = copy.deepcopy(row)
                for col in col_comb:
                    new_row[col] = "UNKNOWN"
                new_rows.append(new_row)
    return pd.DataFrame(new_rows, columns=df.columns)

def count_non_unknown(row):
    """
    Count the number of non-UNKNOWN elements in a given row.

    :param row: The row to be analyzed.
    :return: The number of non-UNKNOWN elements in the row.
    """
    return sum(row != 'UNKNOWN')

def should_exclude_row_2(row):
    """
    Determines whether a row should be excluded based on the values in the row.

    :param row: The row to check.
    :type row: dict
    :return: True if the row should be excluded, False otherwise.
    :rtype: bool
    """
    non_unknown_inds = [index for index, value in row.items() if value != 'UNKNOWN']
    if len(non_unknown_inds) == 2 and set(non_unknown_inds) != {'MARQUES', 'MODELS'}:
        return True
    return False

def should_exclude_row_resolution(row):
    """
    :param row: A dictionary representing a row in a dataset.
    :return: A boolean indicating whether the row should be excluded based on the resolution value.

    The `should_exclude_row_resolution` method checks if the resolution value of a given row is 'UNKNOWN'. If it is, the method returns True, indicating that the row should be excluded.
    * Otherwise, it returns False.
    """
    if row['RESOLUTION'] == 'UNKNOWN':
        return True
    return False

def should_exclude_row_3(row):
    """
    Determines whether a given row should be excluded based on specific conditions.

    :param row: A dictionary representing a row with column names as keys and corresponding values.
    :type row: dict
    :return: A boolean indicating whether the row should be excluded (True) or not (False).
    :rtype: bool
    """
    non_unknown_inds = [index for index, value in row.items() if value != 'UNKNOWN']
    if len(non_unknown_inds) == 2 and set(non_unknown_inds) != {'MARQUES', 'SPECTRUM_TYPE'}:
        return True
    return False

def should_exclude_row_4(row):
    """
    Determines whether the given row should be excluded.

    :param row: The row to be evaluated.
    :type row: dict

    :return: True if the row should be excluded, False otherwise.
    :rtype: bool
    """
    non_unknown_inds = [index for index, value in row.items() if value != 'UNKNOWN']
    if len(non_unknown_inds) == 2 and set(non_unknown_inds) != {'MARQUES', 'INSTRUMENT_TYPE'}:
        return True
    return False

def should_exclude_row_6(row):
    """
    Determines if a given row should be excluded, based on certain conditions.

    :param row: A dictionary representing a row of data.
    :type row: dict
    :return: True if the row should be excluded, False otherwise.
    :rtype: bool
    """
    non_unknown_inds = [index for index, value in row.items() if value != 'UNKNOWN']
    if len(non_unknown_inds) == 2 and set(non_unknown_inds) != {'MARQUES', 'RESOLUTION'}:
        return True
    return False

# def should_exclude_row_5(row):
#     """
#     Determines whether row 5 should be excluded based on its values.
#
#     :param row: A dictionary containing the values of row 5.
#     :type row: dict
#     :return: A boolean indicating whether row 5 should be excluded.
#     :rtype: bool
#     """
#     non_unknown_inds = [index for index, value in row.items() if value != 'UNKNOWN']
#     if len(non_unknown_inds) == 2 and set(non_unknown_inds) != {'MARQUES', 'IONISATION'}:
#         return True
#     return False

def should_exclude_row_7(row):
    """
    Determines if row 7 should be excluded.

    :param row: The row of data to evaluate.
    :type row: dict
    :return: True if row 7 should be excluded, False otherwise.
    :rtype: bool
    """
    if row['MODELS'] == 'UNKNOWN' and row['INSTRUMENT_TYPE'] == 'UNKNOWN':
        return True
    return False

def should_exclude_row_8(row):
    """
    :param row: a dictionary representing a row of data
    :return: a boolean indicating whether the row should be excluded or not

    This method takes in a row of data, represented as a dictionary, and checks whether the values of the 'MODELS' and 'IONISATION' keys are both 'UNKNOWN'. If they are, it returns True
    *, indicating that the row should be excluded. Otherwise, it returns False, indicating that the row should not be excluded.

    Example usage:
    row = {'MODELS': 'UNKNOWN', 'IONISATION': 'UNKNOWN'}
    exclude = should_exclude_row_8(row)
    print(exclude)  # Output: True
    """
    if row['MODELS'] == 'UNKNOWN' and row['IONISATION'] == 'UNKNOWN':
        return True
    return False

# lire toutes les feuilles du fichier Excel
dfs = pd.read_excel(r".\Instrument_tree.xlsx", sheet_name=None)

all_dfs = []  # une liste pour stocker tous les DataFrames filtrés

# dfs est un dictionnaire avec nom de la feuille de calcul en tant que clé et df comme valeur
for sheet_name, df in dfs.items():
    df = df.astype(str)

    # Générer le DataFrame avec les combinaisons manquantes
    df_missing_comb_multi = gen_missing_comb_multi(df)

    df_missing_comb_multi = pd.DataFrame(df_missing_comb_multi, columns=["MARQUES", "MODELS", "SPECTRUM_TYPE", "INSTRUMENT_TYPE", "IONISATION", "RESOLUTION"])

    # Apply this function to each row
    non_unknown_counts = df_missing_comb_multi.apply(count_non_unknown, axis=1)

    # Filter out rows where only 'MODELS' or no columns are not 'UNKNOWN'
    df_missing_comb_multi = df_missing_comb_multi[(non_unknown_counts > 1) | ((non_unknown_counts == 1) & (df_missing_comb_multi['MODELS'] != 'UNKNOWN'))]

    rows_to_exclude = df_missing_comb_multi.apply(should_exclude_row_2, axis=1)
    df_filtered = df_missing_comb_multi.loc[~rows_to_exclude]

    rows_to_exclude_resolution = df_filtered.apply(should_exclude_row_resolution, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude_resolution]

    rows_to_exclude = df_filtered.apply(should_exclude_row_3, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude]

    rows_to_exclude = df_filtered.apply(should_exclude_row_4, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude]

    # rows_to_exclude = df_filtered.apply(should_exclude_row_5, axis=1)
    # df_filtered = df_filtered.loc[~rows_to_exclude]

    rows_to_exclude = df_filtered.apply(should_exclude_row_6, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude]

    rows_to_exclude = df_filtered.apply(should_exclude_row_7, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude]

    rows_to_exclude = df_filtered.apply(should_exclude_row_8, axis=1)
    df_filtered = df_filtered.loc[~rows_to_exclude]

    # ajoute le DataFrame filtré à la liste
    all_dfs.append(df_filtered)

# Concaténez tous les DataFrames en un seul
final_df = pd.concat(all_dfs)

# Écrivez le DataFrame final dans le fichier Excel
final_df.to_excel(r".\permutations.xlsx", index=False)





