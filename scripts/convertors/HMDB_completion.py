import pandas as pd

global HMDB_df
HMDB_df_1 = pd.read_csv("../datas/HMDB_part1.csv", sep=";", encoding="UTF-8")
HMDB_df_2 = pd.read_csv("../datas/HMDB_part2.csv", sep=";", encoding="UTF-8")
HMDB_df_3 = pd.read_csv("../datas/HMDB_part3.csv", sep=";", encoding="UTF-8")
HMDB_df_4 = pd.read_csv("../datas/HMDB_part4.csv", sep=";", encoding="UTF-8")
HMDB_df = pd.concat([HMDB_df_1, HMDB_df_2, HMDB_df_3, HMDB_df_4])

def complete_HMDB_spectrum(spectrum):
    """
    Retrieves additional information from HMDB database for a given spectrum.

    :param spectrum: A dictionary representing the spectrum.
        The dictionary should contain a key named "database-id", which indicates
        the HMDB ID of the spectrum.
    :return: The spectrum dictionary with additional information retrieved from HMDB database.
        The dictionary will be updated with the following keys: "INCHIKEY", "SMILES", "INCHI",
        "MOLECULAR_FORMULA", and "ACC_MASS".
        If the HMDB ID is not found in the HMDB database, the method will return None.
    """
    if "database-id" in spectrum:
        HMDB_ID = spectrum["database-id"]

        line = HMDB_df.loc[HMDB_df['EXTERNAL_ID'] == HMDB_ID]

        if not line.empty:
            spectrum["INCHIKEY"] = line["INCHIKEY"].values[0]
            spectrum["SMILES"] = line["SMILES"].values[0]
            spectrum["INCHI"] = line["INCHI"].values[0]
            spectrum["MOLECULAR_FORMULA"] = line["MOLECULAR_FORMULA"].values[0]
            return spectrum

    return None



