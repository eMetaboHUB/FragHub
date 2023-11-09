from upsetplot import *
from matplotlib import pyplot
import pandas as pd
import os

path= r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\Oct2023"

concatenated_df = pd.DataFrame()

for files in os.listdir(path):
    if files.endswith(".csv") and files != "concats.csv":
        df = pd.read_csv(os.path.join(path,files), sep=";", encoding="UTF-8")

        concatenated_df = pd.concat([concatenated_df, df], axis=0)

concatenated_df = concatenated_df[['FILENAME', 'INCHIKEY']]

# Utilisez la méthode rename() avec des expressions régulières
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*MSMS_.*', 'MSMS_Public', regex=True)
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*MoNA.*', 'MoNA', regex=True)
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*XML_.*', 'HMDB', regex=True)
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*MassBank.*', 'MassBank', regex=True)
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*ALL_GNPS.*', 'GNPS', regex=True)
concatenated_df["FILENAME"] = concatenated_df["FILENAME"].replace(r'.*HMDB.*', 'HMDB', regex=True)

# concatenated_df.to_csv(r'C:\Users\Axel\Documents\Présentations\MSP\datas diagram\Oct2023\concats.csv',sep=";",encoding="UTF-8", index=False)

# Utilisez pivot_table pour réorganiser les données
concatenated_df = concatenated_df.pivot_table(index='INCHIKEY', columns='FILENAME', aggfunc='first', fill_value='')


# Supprimez les doublons
df = concatenated_df.drop_duplicates(subset=['FILENAME', 'INCHIKEY'])

# Utilisez la fonction pivot pour réorganiser les données
pivot_df = df.pivot(index='INCHIKEY', columns='FILENAME', values='INCHIKEY').reset_index()

# Remplacez les valeurs non définies par une chaîne vide
pivot_df = pivot_df.fillna('')

pivot_df = pivot_df.drop('INCHIKEY', axis=1)

print(pivot_df)

UpSet(pivot_df)


