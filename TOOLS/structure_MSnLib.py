import os
import pandas as pd
import json

# Chemin du dossier contenant les fichiers .txt
dossier = r"D:\Axel\DATABASE\MSnLib"

# Lister tous les fichiers dans le dossier
fichiers = os.listdir(dossier)

# Filtrer les fichiers .txt
fichiers_txt = [f for f in fichiers if f.endswith('.txt')]

# Liste pour stocker tous les DataFrames
dataframes = []

# Lire chaque fichier .txt avec pandas et ajouter à la liste des DataFrames
for fichier in fichiers_txt:
    chemin_fichier = os.path.join(dossier, fichier)
    df = pd.read_csv(chemin_fichier, sep=';')
    dataframes.append(df)

# Concaténer tous les DataFrames
df_concat = pd.concat(dataframes, ignore_index=True)

# Renommer la colonne 'MS2Peaks' en 'peaks'
df_concat.rename(columns={'MS2Peaks': 'peaks_json'}, inplace=True)


# Transformer chaque valeur de 'peaks' en matrice
def transform_peaks(peaks_str):
    peak_list = []
    for peak in peaks_str.strip().split('\n'):
        mass, intensity = map(float, peak.split('\t'))
        peak_list.append([mass, intensity])
    return peak_list


df_concat['peaks_json'] = df_concat['peaks_json'].apply(transform_peaks)

# Convertir en dictionnaire format records
dict_records = df_concat.to_dict(orient='records')

# Sauvegarder en fichier JSON
chemin_sortie_json = os.path.join(dossier, "output.json")
with open(chemin_sortie_json, 'w', encoding='utf-8') as f:
    json.dump(dict_records, f, ensure_ascii=False, indent=4)

print(f"Le fichier JSON a été enregistré dans : {chemin_sortie_json}")
