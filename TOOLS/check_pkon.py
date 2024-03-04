import pandas as pd
import os

path = r"C:/Users/Axel/Downloads"

file = os.path.join(path,"report(1).csv")

df = pd.read_csv(file,sep=",",quotechar='"', encoding="latin1")

df_filtré = pd.DataFrame()

# Parcourez chaque ligne du DataFrame d'origine (df)
for index, row in df.iterrows():
    # Parcourez chaque colonne (cellule) dans la ligne
    for column in df.columns:
        # Vérifiez si la cellule a plus de 4000 caractères
        if len(str(row[column])) > 3998:
            # Si oui, ajoutez la ligne entière au DataFrame filtré
            df_filtré = pd.concat([df_filtré, row.to_frame().T], ignore_index=True)

print(df_filtré)
