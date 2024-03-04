import pandas as pd
import re


def csv_clean_to_msp(CSV_df):
    POS_LC = []
    for index, row in CSV_df.iterrows():
        SPECTRUM = ""
        SPECTRUM = SPECTRUM + "FILENAME: " + row["FILENAME"] + "\n"
        SPECTRUM = SPECTRUM + "PREDICTED: " + row["PREDICTED"] + "\n"
        SPECTRUM = SPECTRUM + "SPECTRUMID: " + row["SPECTRUMID"] + "\n"
        SPECTRUM = SPECTRUM + "RESOLUTION: " + row["RESOLUTION"] + "\n"
        SPECTRUM = SPECTRUM + "SYNON: " + row["SYNON"] + "\n"
        SPECTRUM = SPECTRUM + "CHARGE: " + row["CHARGE"] + "\n"
        SPECTRUM = SPECTRUM + "PARENTMASS: " + row["PARENTMASS"] + "\n"
        SPECTRUM = SPECTRUM + "IONIZATION: " + row["IONIZATION"] + "\n"
        SPECTRUM = SPECTRUM + "MSLEVEL: " + row["MSLEVEL"] + "\n"
        SPECTRUM = SPECTRUM + "FRAGMENTATIONMODE: " + row["FRAGMENTATIONMODE"] + "\n"
        SPECTRUM = SPECTRUM + "NAME: " + row["NAME"] + "\n"
        SPECTRUM = SPECTRUM + "PRECURSORMZ: " + row["PRECURSORMZ"] + "\n"
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


    content = "\n\n".join(POS_LC)

    with open(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\COSIN\TEST_colision_temp.msp","w",encoding="UTF-8") as buffer:
        buffer.write(content)

CSV =  r"C:\Users\Axel\Documents\KNIME\MSP3_check\CSV\FINAL_POS\POS_LC_clean.csv"

CSV_df = pd.read_csv(CSV, sep=";", encoding="UTF-8")

# Convertir toutes les colonnes en chaînes de caractères
CSV_df = CSV_df.astype(str)

# PETIT NETTOYAGE (juste pour les heatmaps)
CSV_df = CSV_df[CSV_df['COLLISIONENERGY'] != "None"]
CSV_df = CSV_df[CSV_df['COLLISIONENERGY'] != "nan"]
CSV_df = CSV_df[CSV_df['INSTRUMENT'] != "None"]
CSV_df = CSV_df[CSV_df['INSTRUMENT'] != "N/A"]


csv_clean_to_msp(CSV_df)