import pandas as pd
import os

NEG_path = r"C:\Users\Axel\Documents\KNIME\MSP3_check\CSV\FINAL_NEG"
POS_path = r"C:\Users\Axel\Documents\KNIME\MSP3_check\CSV\FINAL_POS"


# List of columns to keep
columns_to_keep = ["INCHI", "INCHIKEY", "SMILES"]

# POS
for files in os.listdir(POS_path):
    if files.endswith(".csv"):
        try:
            df_temp = pd.read_csv(os.path.join(POS_path,files),sep=";",encoding="UTF-8", usecols=columns_to_keep)
            # Concatenate the current DataFrame with DF_FINAL
            df_temp = df_temp.drop_duplicates(subset=["INCHI", "INCHIKEY", "SMILES"], keep="first")
            df_temp.to_csv(os.path.join(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean",files), sep=";", index=False, encoding="UTF-8")
        except:
            continue

# NEG
for files in os.listdir(NEG_path):
    if files.endswith(".csv"):
        try:
            df_temp = pd.read_csv(os.path.join(NEG_path,files),sep=";",encoding="UTF-8", usecols=columns_to_keep)
            # Concatenate the current DataFrame with DF_FINAL
            df_temp = df_temp.drop_duplicates(subset=["INCHI", "INCHIKEY", "SMILES"], keep="first")
            df_temp.to_csv(os.path.join(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean",files), sep=";", index=False, encoding="UTF-8")
        except:
            continue




