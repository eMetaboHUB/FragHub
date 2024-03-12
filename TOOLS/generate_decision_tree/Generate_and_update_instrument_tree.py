from collections import defaultdict
import pandas as pd
import itertools
import copy
import json

def gen_missing_comb_multi(df):
    """
    Generate missing permutations for multiple columns.

    :param df: The DataFrame containing the columns.
    :return: A new DataFrame with additional rows representing the missing permutations.

    """
    new_rows = []

    # Get a list of all the column names
    columns = df.columns.tolist()

    for r in range(len(columns) + 1):
        # Generate all combinations of a given length
        for subset in itertools.combinations(columns, r):
            for _, row in df.iterrows():
                new_row = row.copy()
                # Set the selected columns to "UNKNOWN"
                for col in subset:
                    new_row[col] = "UNKNOWN"
                new_rows.append(new_row.tolist())

    return pd.DataFrame(new_rows, columns=df.columns)

def format_solution(row, non_unknown_marques):
    """
    :param row: dictionary containing various attributes
    :param non_unknown_marques: string representing non-unknown marques
    :return: formatted solution based on given conditions

    """
    # complete la marque si on a la modele
    MARQUES = non_unknown_marques if row['MARQUES'] == 'UNKNOWN' and row['MODELS'] != 'UNKNOWN' else row['MARQUES']

    # complete 'marque instrument' si on a pas de nom de model
    MODELS = 'instrument' if row['MODELS'] == 'UNKNOWN' and row['MARQUES'] != 'UNKNOWN' else row['MODELS']

    # complete LC ou GC si on à la type d'ionisation
    if row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] != 'UNKNOWN':
        if row['IONISATION'] in ['ESI', 'APPI', 'APCI', 'MALDI', 'FAB', 'FD', 'EI']:
            SPECTRUM_TYPE = 'LC'
        elif row['IONISATION'] in ['CI', 'PTR', 'EI']:
            SPECTRUM_TYPE = 'GC'
        else:
            SPECTRUM_TYPE = row['SPECTRUM_TYPE']
    else:
        SPECTRUM_TYPE = row['SPECTRUM_TYPE']

    if row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row["IONISATION"]}-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {row['SPECTRUM_TYPE']}-UNKNOWN-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, {SPECTRUM_TYPE}-{row["IONISATION"]}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {MODELS}, UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, {SPECTRUM_TYPE}-{row["IONISATION"]}-{row['INSTRUMENT_TYPE_2']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{MARQUES} {row['MODELS']}, UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"UNKNOWN, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"UNKNOWN, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"UNKNOWN, {row['SPECTRUM_TYPE']}-{row["IONISATION"]}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"UNKNOWN, {row['SPECTRUM_TYPE']}-UNKNOWN-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"UNKNOWN, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"UNKNOWN, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"UNKNOWN, {SPECTRUM_TYPE}-{row["IONISATION"]}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] == 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"UNKNOWN, UNKNOWN, UNKNOWN"

def get_first_non_unknown_marques(df):
    """
    Get the first non-UNKNOWN marque from a DataFrame.

    :param df: The DataFrame containing the 'MARQUES' column.
    :type df: pandas.DataFrame
    :return: The first non-UNKNOWN marque if found, otherwise 'UNKNOWN'.
    :rtype: str
    """
    for marque in df['MARQUES']:
        if marque != 'UNKNOWN':
            return marque
    # If there are no non-UNKNOWN marques, just return 'UNKNOWN'
    return 'UNKNOWN'

def remove_selected_unknowns(df):
    """
    Remove all rows that have 'UNKNOWN' in the selected columns.

    :param df: The DataFrame containing the columns.
    :return: A new DataFrame without rows with 'UNKNOWN' in the selected columns.
    """
    cols_to_check = ['MARQUES', 'MODELS', 'SPECTRUM_TYPE', 'INSTRUMENT_TYPE', 'IONISATION']
    return df[~(df[cols_to_check] == 'UNKNOWN').all(axis=1)]

def exclude_useless_rows(df):
    """
    Remove useless rows from the given DataFrame.

    :param df: The DataFrame to remove useless rows from.
    :type df: pandas.DataFrame
    :return: The DataFrame with useless rows removed.
    :rtype: pandas.DataFrame
    """
    df = remove_selected_unknowns(df)

    return df

def process_keys(obj):
    """
    Process the keys within a given object.

    :param obj: The object to be processed.
    :return: The object with the processed keys.
    """
    if isinstance(obj, dict):
        return {key if key == 'SOLUTION' else key.lower().strip(): process_keys(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [process_keys(item) for item in obj]
    else:
        return obj

def tree():
    """
    Return a `defaultdict` object initialized to a nested tree structure.

    :return: A `defaultdict` object with each key mapping to another nested `defaultdict` object.
    """
    return defaultdict(tree)

def dictify(nested_dict):
    """
    :param nested_dict: A nested dictionary containing key-value pairs.
    :return: A flattened dictionary where all the keys from the nested dictionary are merged into one, and the values remain the same.

    """
    result = {}
    for key, value in nested_dict.items():
        if key == "SOLUTION":
            result[key] = value
        elif isinstance(value, dict):
            result[key] = dictify(value)
    return result

def fill_unknown_instruments_type(group):
    """

    :param group: A pandas Series or DataFrame column containing instrument types. It can have the value 'UNKNOWN' for unknown instrument types.
    :return: The input group after replacing all occurrences of 'UNKNOWN' with the most common non-'UNKNOWN' instrument type.

    """
    mask = group == 'UNKNOWN'
    if not group[~mask].empty:
        # Find the mode(s) of the group
        modes = group[~mask].mode()
        # Convert modes to strings and sort them by length
        selected_value = modes.astype(str).sort_values(key=lambda x: x.str.len()).iloc[0]
        group[mask] = selected_value
    return group

if __name__ == '__main__':
    # lire toutes les feuilles du fichier Excel
    dfs = pd.read_excel(r".\Instrument_tree.xlsx", sheet_name=None)

    all_dfs = []  # une liste pour stocker tous les DataFrames filtrés

    # dfs est un dictionnaire avec nom de la feuille de calcul en tant que clé et df comme valeur
    for sheet_name, df in dfs.items():
        df = df.astype(str)

        # Appliquer strip à toutes les valeurs de DataFrame
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        # Générer le DataFrame avec les combinaisons manquantes
        df_missing_comb_multi = gen_missing_comb_multi(df)

        df_missing_comb_multi = pd.DataFrame(df_missing_comb_multi, columns=["MARQUES", "MODELS", "SPECTRUM_TYPE", "INSTRUMENT_TYPE", "IONISATION", "RESOLUTION"])

        df_filtered = exclude_useless_rows(df_missing_comb_multi).copy()


        df_filtered['INSTRUMENT_TYPE_2'] = df_filtered.groupby('MODELS')['INSTRUMENT_TYPE'].transform(fill_unknown_instruments_type)

        non_unknown_marques = get_first_non_unknown_marques(df_filtered)

        df_filtered['SOLUTION'] = df_filtered.apply(lambda row: format_solution(row, non_unknown_marques), axis=1)

        df_filtered.drop('INSTRUMENT_TYPE_2', axis=1, inplace=True)

        # Modifier certaines valeurs
        # Modifier les valeurs dans la colonne MARQUES
        df_filtered['MARQUES'] = df_filtered['MARQUES'].replace({'AB SCIEX': 'SCIEX', 'Thermo Fisher Scientific': 'Thermo'})

        df_filtered = df_filtered.apply(lambda x: x.str.lower() if (x.dtype == "object" and x.name != 'SOLUTION') else x)

        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("-tof", "tof", regex=True)
        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("q-", "q", regex=True)
        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("q exactive", " qexactive ", regex=True)
        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("-", " ", regex=True)
        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("triple(-| )?tof", " qqq ", regex=True)
        df_filtered['MODELS'] = df_filtered['MODELS'].str.replace("triple(-| )?quad", " qqq ", regex=True)

        # Appliquer strip à toutes les valeurs de DataFrame
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        # trie le DataFrame en fonction de la longueur de la chaîne 'MODELS'
        df_filtered = df_filtered.sort_values(by='MODELS', key=lambda x: x.str.len(), ascending=False)

        # ajoute le DataFrame filtré à la liste
        all_dfs.append(df_filtered)

        # Concaténez tous les DataFrames en un seul
        final_df = pd.concat(all_dfs)

    # Écrivez le DataFrame final dans le fichier Excel
    final_df.to_excel(r".\GRAPH_CONSTRUCT_SOLUTIONS.xlsx", index=False)

    root = tree()
    for _, row in final_df.iterrows():
        node = root
        for val in row[:-1]:  # Ignore the last value, which is "SOLUTION"
            node = node[val]
        node['SOLUTION'] = row['SOLUTION']

    D = dictify(root)
    D = process_keys(D)  # New line to process keys.

    with open(r'.\instruments_tree.json', 'w') as f:
        json.dump(D, f, indent=2)







