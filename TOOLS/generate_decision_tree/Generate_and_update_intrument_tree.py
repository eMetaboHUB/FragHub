import pandas as pd
import itertools
import copy

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
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] != 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row["IONISATION"]}-UNKNOWN, {row['RESOLUTION']}"
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
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] != 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {row['SPECTRUM_TYPE']}-UNKNOWN-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] != 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, UNKNOWN-UNKNOWN-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] != 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, {SPECTRUM_TYPE}-{row["IONISATION"]}-UNKNOWN, {row['RESOLUTION']}"
    elif row["MARQUES"] == 'UNKNOWN' and row["MODELS"] != 'UNKNOWN' and row["SPECTRUM_TYPE"] == 'UNKNOWN' and row["IONISATION"] == 'UNKNOWN' and row["INSTRUMENT_TYPE"] == 'UNKNOWN':
        return f"{row["MARQUES"]} {row['MODELS']}, UNKNOWN, {row['RESOLUTION']}"
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

if __name__ == '__main__':
    # lire toutes les feuilles du fichier Excel
    dfs = pd.read_excel(r".\Instrument_tree.xlsx", sheet_name=None)

    all_dfs = []  # une liste pour stocker tous les DataFrames filtrés

    # dfs est un dictionnaire avec nom de la feuille de calcul en tant que clé et df comme valeur
    for sheet_name, df in dfs.items():
        df = df.astype(str)

        # Générer le DataFrame avec les combinaisons manquantes
        df_missing_comb_multi = gen_missing_comb_multi(df)

        df_missing_comb_multi = pd.DataFrame(df_missing_comb_multi, columns=["MARQUES", "MODELS", "SPECTRUM_TYPE", "INSTRUMENT_TYPE", "IONISATION", "RESOLUTION"])

        df_filtered = exclude_useless_rows(df_missing_comb_multi).copy()

        non_unknown_marques = get_first_non_unknown_marques(df_filtered)

        # ici competer les infos manquantes


        df_filtered['SOLUTION'] = df_filtered.apply(lambda row: format_solution(row, non_unknown_marques), axis=1)

        # ajoute le DataFrame filtré à la liste
        all_dfs.append(df_filtered)

        # Concaténez tous les DataFrames en un seul
        final_df = pd.concat(all_dfs)

    # Écrivez le DataFrame final dans le fichier Excel
    final_df.to_excel(r".\GRAPH_CONSTRUCT_SOLUTIONS.xlsx", index=False)





