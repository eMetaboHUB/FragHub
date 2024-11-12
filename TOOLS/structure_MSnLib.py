import os
import pandas as pd
import json
from tqdm import tqdm

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
    print(f"reading: {chemin_fichier}")
    df = pd.read_csv(chemin_fichier, sep=';')
    dataframes.append(df)

# Concaténer tous les DataFrames
df_concat = pd.concat(dataframes, ignore_index=True)

# Renommer la colonne 'MS2Peaks' en 'peaks_json'
df_concat.rename(columns={'MS2Peaks': 'peaks_json'}, inplace=True)

# Remplacer les cellules vides par des chaînes vides
df_concat.fillna("", inplace=True)


# Transformer chaque valeur de 'peaks_json' en matrice
def transform_peaks(peaks_str):
    peak_list = []
    for peak in str(peaks_str).strip().split('\n'):
        parts = peak.split('\t')
        if len(parts) == 2:
            try:
                mass, intensity = map(float, parts)
                peak_list.append([mass, intensity])
            except ValueError:
                print(f"Skipping invalid peak: {peak}")
        elif len(parts) == 1 and parts[0]:  # Gérer le cas où il n'y a qu'une seule valeur non vide
            try:
                mass = float(parts[0])
                peak_list.append([mass, 0])
            except ValueError:
                print(f"Skipping invalid peak: {peak}")
        else:
            print(f"Skipping invalid or empty peak: {peak}")
    return peak_list


# Utilisation de la barre de progression avec tqdm pour la méthode apply
tqdm.pandas()

# Application de la transformation avec barre de progression
df_concat['peaks_json'] = df_concat['peaks_json'].progress_apply(transform_peaks)

# Convertir en dictionnaire format records
print("df to dict")
dict_records = df_concat.to_dict(orient='records')

# Sauvegarder en fichier JSON avec une barre de progression
chemin_sortie_json = os.path.join(dossier, "MetaboAnalyste.json")
print("Écriture des fichiers JSON")
with open(chemin_sortie_json, 'w', encoding='utf-8') as f:
    f.write('[\n')  # Début de la liste JSON
    # Écriture avec barre de progression
    for i, record in enumerate(tqdm(dict_records, desc="Sauvegarde en JSON")):
        f.write(json.dumps(record, ensure_ascii=False))
        if i < len(dict_records) - 1:
            f.write(',\n')  # Ajoute une virgule après chaque élément sauf le dernier
    f.write('\n]')  # Fin de la liste JSON

print(f"Le fichier JSON a été enregistré dans : {chemin_sortie_json}")
