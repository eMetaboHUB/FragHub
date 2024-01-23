import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Étape 1 : Lire les données CSV dans un DataFrame
df_1 = pd.read_csv(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean\NEG.csv", sep=";", encoding="UTF-8")
df_2 = pd.read_csv(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean\POS.csv", sep=";", encoding="UTF-8")

df = pd.concat([df_1, df_2], ignore_index=True)
df = df.sort_values(by='FILENAME')

df = df.dropna(subset=["INCHIKEY"])

# df.to_csv(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean\FULL.csv",sep=";",encoding="UTF-8",index=False)

# Étape 2 : Préparez les données pour le graphique à barres doubles
filename = df["FILENAME"].unique()
inchikey = df["INCHIKEY"]

# Création d'un index personnalisé pour chaque DB
index = np.arange(len(filename))

# Calcul du nombre total d'inchikeys pour chaque DB
n_inchikey_tot = df.groupby("FILENAME")["INCHIKEY"].count()

# Calcul du nombre d'inchikeys uniques pour chaque DB
n_inchikey_uni = df.groupby("FILENAME")["INCHIKEY"].nunique()

# Largeur des barres
bar_width = 0.4

# Étape 3 : Créez le graphique à barres doubles avec les barres séparées
plt.figure(figsize=(16, 10))  # Facultatif : définir la taille du graphique

# Création de la première barre "Nombre total d'inchikeys" (au-dessus de la deuxième)
bars1 = plt.barh(index + bar_width/2, n_inchikey_tot, bar_width, label="Total number of inchikeys")
bars2 = plt.barh(index - bar_width/2, n_inchikey_uni, bar_width, label="Total number of unique inchikeys", color='orange')

# Définir les étiquettes des axes
plt.ylabel("DB")
plt.xlabel("Number of inchikeys")

# Définir le titre du graphique

# Définir les étiquettes de l'axe y
plt.yticks(index, filename)

# Ajouter les nombres d'inchikeys au-dessus des barres
for bar in bars1 + bars2:
    width = bar.get_width()
    yloc = bar.get_y() + bar.get_height() / 2
    plt.text(width, yloc, f'{int(width)}', ha='left', va='center', color='black')

# Ajuster la mise en page
plt.tight_layout()

# Afficher le graphique
plt.legend()
plt.show()
