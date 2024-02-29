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
    Format the solution based on the given row and non_unknown_marques.

    :param row: A dictionary representing a row of data.
    :type row: dict
    :param non_unknown_marques: A string representing the non-UNKNOWN marques.
    :type non_unknown_marques: str
    :return: The formatted solution.
    :rtype: str
    """
    marque = non_unknown_marques if row['MARQUES'] == 'UNKNOWN' and row['MODELS'] != 'UNKNOWN' else row['MARQUES']
    model = 'instrument' if row['MODELS'] == 'UNKNOWN' and row['MARQUES'] != 'UNKNOWN' else row['MODELS']

    if row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] != 'UNKNOWN':
        if row['IONISATION'] in ['ESI', 'APPI', 'APCI', 'MALDI', 'FAB', 'FD', 'EI']:
            spectrum_type = 'LC'
        elif row['IONISATION'] in ['CI', 'PTR', 'EI']:
            spectrum_type = 'GC'
        else:
            spectrum_type = row['SPECTRUM_TYPE']
    else:
        spectrum_type = row['SPECTRUM_TYPE']

    if row['SPECTRUM_TYPE'] != 'UNKNOWN' and row['IONISATION'] != 'UNKNOWN' and row['INSTRUMENT_TYPE'] == 'UNKNOWN':
        return f"{marque} {model}, {spectrum_type}-{row['IONISATION']}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] == 'UNKNOWN' and row['INSTRUMENT_TYPE'] != 'UNKNOWN':
        return f"UNKNOWN, {row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] != 'UNKNOWN' and row['IONISATION'] != 'UNKNOWN' and row['INSTRUMENT_TYPE'] != 'UNKNOWN':
        return f"{marque} {model}, {row['SPECTRUM_TYPE']}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] == 'UNKNOWN' and row['INSTRUMENT_TYPE'] != 'UNKNOWN':
        return f"{marque} {model}, {spectrum_type}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] != 'UNKNOWN' and row['IONISATION'] == 'UNKNOWN' and row['INSTRUMENT_TYPE'] == 'UNKNOWN':
        return f"{marque} {model}, {spectrum_type}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] != 'UNKNOWN' and row['INSTRUMENT_TYPE'] == 'UNKNOWN':
        return f"{marque} {model}, {spectrum_type}-{row['IONISATION']}, {row['RESOLUTION']}"
    elif row['SPECTRUM_TYPE'] == 'UNKNOWN' and row['IONISATION'] == 'UNKNOWN' and row['INSTRUMENT_TYPE'] == 'UNKNOWN':
        return f"{marque} {model}, UNKNOWN, {row['RESOLUTION']}"
    elif marque == 'UNKNOWN' and model == 'UNKNOWN':
        return f"UNKNOWN, {spectrum_type}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"
    else:
        return f"{marque} {model}, {spectrum_type}-{row['IONISATION']}-{row['INSTRUMENT_TYPE']}, {row['RESOLUTION']}"

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

if __name__ == '__main__':
    # lire toutes les feuilles du fichier Excel
    dfs = pd.read_excel(r".\TEST.xlsx", sheet_name=None)

    all_dfs = []  # une liste pour stocker tous les DataFrames filtrés

    # dfs est un dictionnaire avec nom de la feuille de calcul en tant que clé et df comme valeur
    for sheet_name, df in dfs.items():
        df = df.astype(str)

        # Générer le DataFrame avec les combinaisons manquantes
        df_missing_comb_multi = gen_missing_comb_multi(df)

        df_missing_comb_multi = pd.DataFrame(df_missing_comb_multi, columns=["MARQUES", "MODELS", "SPECTRUM_TYPE", "INSTRUMENT_TYPE", "IONISATION", "RESOLUTION"])

        non_unknown_marques = get_first_non_unknown_marques(df_missing_comb_multi)

        # rajoutons 'SOLUTION' à chaque ligne
        non_unknown_marques = get_first_non_unknown_marques(df_missing_comb_multi)
        df_missing_comb_multi['SOLUTION'] = df_missing_comb_multi.apply(lambda row: format_solution(row, non_unknown_marques), axis=1)

        # ajoute le DataFrame filtré à la liste
        all_dfs.append(df_missing_comb_multi)

        # Concaténez tous les DataFrames en un seul
        final_df = pd.concat(all_dfs)

    # Écrivez le DataFrame final dans le fichier Excel
    final_df.to_excel(r".\GRAPH_CONSTRUCT_SOLUTIONS_TEST.xlsx", index=False)





