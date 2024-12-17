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
    return f'FILENAME={DF_row["FILENAME"] if DF_row["FILENAME"] else "UNKNOWN"}; PREDICTED={DF_row["PREDICTED"] if DF_row["PREDICTED"] else "UNKNOWN"}; SPLASH={DF_row["SPLASH"] if DF_row["SPLASH"] else "UNKNOWN"}; SPECTRUMID={DF_row["SPECTRUMID"] if DF_row["SPECTRUMID"] else "UNKNOWN"}; RESOLUTION={DF_row["RESOLUTION"] if DF_row["RESOLUTION"] else "UNKNOWN"}; SYNON={DF_row["SYNON"] if DF_row["SYNON"] else "UNKNOWN"}; FRAGMENTATIONMODE={DF_row["FRAGMENTATIONMODE"] if DF_row["FRAGMENTATIONMODE"] else "UNKNOWN"}; AVERAGEMASS={DF_row["AVERAGEMASS"] if DF_row["AVERAGEMASS"] else "UNKNOWN"}; ENTROPY={DF_row["ENTROPY"] if DF_row["ENTROPY"] else "UNKNOWN"}; ONTOLOGIES = "CLASSYFIRE_SUPERCLASS={DF_row['CLASSYFIRE_SUPERCLASS'] if DF_row['CLASSYFIRE_SUPERCLASS'] else 'UNKNOWN'}, CLASSYFIRE_CLASS = {DF_row['CLASSYFIRE_CLASS'] if DF_row['CLASSYFIRE_CLASS'] else 'UNKNOWN'}, CLASSYFIRE_SUBCLASS = {DF_row['CLASSYFIRE_SUBCLASS'] if DF_row['CLASSYFIRE_SUBCLASS'] else 'UNKNOWN'}, NPCLASS_PATHWAY = {DF_row['NPCLASS_PATHWAY'] if DF_row['NPCLASS_PATHWAY'] else 'UNKNOWN'}, NPCLASS_SUPERCLASS = {DF_row['NPCLASS_SUPERCLASS'] if DF_row['NPCLASS_SUPERCLASS'] else 'UNKNOWN'}, NPCLASS_CLASS = {DF_row['NPCLASS_CLASS'] if DF_row['NPCLASS_CLASS'] else 'UNKNOWN'}"'

def dataframe_to_msp(dataframe, name, progress_callback=None, total_items_callback=None, prefix_callback=None,
                     item_type_callback=None):
    """
    Convertit un DataFrame contenant des données spectrales en une liste de chaînes formatées représentant les spectres.

    :param dataframe: Un DataFrame pandas contenant les données spectrales.
    :param name: Le nom du spectre (sert de contexte pour les callbacks).
    :param progress_callback: Fonction callable pour notifier la progression.
    :param total_items_callback: Fonction callable pour définir le total des éléments.
    :param prefix_callback: Fonction callable pour signaler le contexte de la tâche.
    :param item_type_callback: Fonction callable pour signaler le type des éléments traités.
    :return: Une liste de chaînes formatées représentant les spectres.
    """

    # Convertir toutes les colonnes en type chaîne
    dataframe = dataframe.astype(str)

    # Indique la tâche en cours
    if prefix_callback:
        prefix_callback(f"Formatting {name} to MSP")

    # Définir le type d'éléments traités
    if item_type_callback:
        item_type_callback("rows")

    # Signaler le total des éléments à traiter
    total_rows = len(dataframe)
    if total_items_callback:
        total_items_callback(total_rows, 0)  # Initialisation des éléments complétés

    # Initialiser une liste pour les chaînes formatées
    spectrum_list = []

    # Parcourir chaque ligne du DataFrame
    for index, row in dataframe.iterrows():
        # Formater le commentaire pour la ligne courante
        COMMENTS = format_comments(row)

        # Créer la chaîne du spectre basé sur les valeurs de la ligne
        SPECTRUM = ""
        SPECTRUM += "NAME: " + (row["NAME"] if row["NAME"] else "UNKNOWN") + "\n"
        SPECTRUM += "PRECURSORMZ: " + (row["PRECURSORMZ"] if row["PRECURSORMZ"] else "UNKNOWN") + "\n"
        SPECTRUM += "PRECURSORTYPE: " + (row["PRECURSORTYPE"] if row["PRECURSORTYPE"] else "UNKNOWN") + "\n"
        SPECTRUM += "FORMULA: " + (row["FORMULA"] if row["FORMULA"] else "UNKNOWN") + "\n"
        SPECTRUM += "INCHIKEY: " + (row["INCHIKEY"] if row["INCHIKEY"] else "UNKNOWN") + "\n"
        SPECTRUM += "INCHI: " + (row["INCHI"] if row["INCHI"] else "UNKNOWN") + "\n"
        SPECTRUM += "SMILES: " + (row["SMILES"] if row["SMILES"] else "UNKNOWN") + "\n"
        SPECTRUM += "RETENTIONTIME: " + (row["RETENTIONTIME"] if row["RETENTIONTIME"] else "UNKNOWN") + "\n"
        SPECTRUM += "IONMODE: " + (row["IONMODE"] if row["IONMODE"] else "UNKNOWN") + "\n"
        SPECTRUM += "INSTRUMENTTYPE: " + (row["INSTRUMENTTYPE"] if row["INSTRUMENTTYPE"] else "UNKNOWN") + "\n"
        SPECTRUM += "INSTRUMENT: " + (row["INSTRUMENT"] if row["INSTRUMENT"] else "UNKNOWN") + "\n"
        SPECTRUM += "COLLISIONENERGY: " + (row["COLLISIONENERGY"] if row["COLLISIONENERGY"] else "UNKNOWN") + "\n"
        SPECTRUM += "EXACTMASS: " + (row["EXACTMASS"] if row["EXACTMASS"] else "UNKNOWN") + "\n"
        SPECTRUM += "IONIZATION: " + (row["IONIZATION"] if row["IONIZATION"] else "UNKNOWN") + "\n"
        SPECTRUM += "MSLEVEL: " + (row["MSLEVEL"] if row["MSLEVEL"] else "UNKNOWN") + "\n"
        SPECTRUM += "COMMENT: " + (COMMENTS if COMMENTS else "UNKNOWN") + "\n"
        SPECTRUM += "NUM PEAKS: " + (row["NUM PEAKS"] if row["NUM PEAKS"] else "UNKNOWN") + "\n"
        SPECTRUM += row["PEAKS_LIST"] + "\n"

        # Ajouter le spectre formaté à la liste des spectres
        spectrum_list.append(SPECTRUM)

        # Mise à jour de la progression
        if progress_callback:
            progress_callback(index + 1)

    # Retourner la liste des spectres formatés
    return spectrum_list


def csv_to_msp(POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
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
    POS_LC = dataframe_to_msp(POS_LC_df, "POS_LC", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Create MSP file from POS_LC insilico dataframe
    POS_LC_insilico = dataframe_to_msp(POS_LC_df_insilico, "POS_LC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Create MSP file from POS_GC dataframe
    POS_GC = dataframe_to_msp(POS_GC_df, "POS_GC")
    # Create MSP file from POS_GC insilico dataframe
    POS_GC_insilico = dataframe_to_msp(POS_GC_df_insilico, "POS_GC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Create MSP file from NEG_LC dataframe
    NEG_LC = dataframe_to_msp(NEG_LC_df, "NEG_LC")
    # Create MSP file from NEG_LC insilico dataframe
    NEG_LC_insilico = dataframe_to_msp(NEG_LC_df_insilico, "NEG_LC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
    # Create MSP file from NEG_GC dataframe
    NEG_GC = dataframe_to_msp(NEG_GC_df, "NEG_GC")
    # Create MSP file from NEG_GC insilico dataframe
    NEG_GC_insilico = dataframe_to_msp(NEG_GC_df_insilico, "NEG_GC_insilico", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Return original dataframes as well as the generated MSP files.
    return POS_LC_df, POS_LC, POS_LC_df_insilico, POS_LC_insilico, POS_GC_df, POS_GC, POS_GC_df_insilico, POS_GC_insilico, NEG_LC_df, NEG_LC, NEG_LC_df_insilico, NEG_LC_insilico, NEG_GC_df, NEG_GC, NEG_GC_df_insilico, NEG_GC_insilico
