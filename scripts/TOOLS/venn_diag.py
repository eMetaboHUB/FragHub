import pandas as pd
from venn import venn
from venn import pseudovenn
import matplotlib.pyplot as plt


# This example script simply outputs the node's input table.

# Étape 1 : Lire les données CSV dans un DataFrame
nom_du_fichier = r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\clean\FULL.csv"
df = pd.read_csv(nom_du_fichier, sep=";", encoding="UTF-8")
# df = df.dropna()
# print(df)

# Remplacer 'FILENAME' par 'In-Silico' lorsque 'PREDICTED' est à True
df.loc[df['PREDICTED'] == " true", 'FILENAME'] = 'In-Silico'

# df.to_csv(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\not_clean\ALL_DB_test.csv",sep=";", encoding="UTF-8")
occurence_des_BD = {}


for i in df.index:
    if df["FILENAME"][i] not in occurence_des_BD:
        occurence_des_BD[ df["FILENAME"][i]] = []
        occurence_des_BD[df["FILENAME"][i]].append(df["INCHIKEY"][i].replace(" ",""))
    else:
        occurence_des_BD[df["FILENAME"][i]].append(df["INCHIKEY"][i].replace(" ",""))

for keys,values in occurence_des_BD.items():
    occurence_des_BD[keys] = set(occurence_des_BD[keys])


# Plot the histogram
#fig = plt.figure()
#pseudovenn(occurence_des_BD, cmap="plasma")

# Assign the figure to the output_view variable
fig = plt.figure()
# venn_diagram = venn(occurence_des_BD,fmt="{percentage:.1f}%", legend_loc="upper left")
venn_diagram = pseudovenn(occurence_des_BD, fmt="{percentage:.1f}%", cmap="plasma", hint_hidden=False, legend_loc="upper left")
# Affichez le diagramme de Venn
plt.show()

# import pandas as pd
# from matplotlib_venn import venn2
# import matplotlib.pyplot as plt
#
# # Étape 1 : Lire les données CSV dans un DataFrame
# nom_du_fichier = r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\not_clean\ALL_DB.csv"
# df = pd.read_csv(nom_du_fichier, sep=";", encoding="UTF-8")
# df = df.dropna()
#
# occurence_des_BD = {}
#
# for i in df.index:
#     if df["filename"][i] not in occurence_des_BD:
#         occurence_des_BD[df["filename"][i]] = set()
#     occurence_des_BD[df["filename"][i]].add(df["inchikey"][i])
#
# # Créez le diagramme de Venn
# plt.figure(figsize=(8, 6))  # Définissez la taille de la figure (facultatif)
# venn_diagram = venn2([occurence_des_BD[key] for key in occurence_des_BD],
#                      set_labels=occurence_des_BD.keys())
#
# # Ajoutez le nombre d'éléments à l'intérieur de chaque zone du diagramme
# for label, subset in venn_diagram.set_labels.items():
#     subset.set_text(len(occurence_des_BD[label]))
# for subset in venn_diagram.subset_labels:
#     subset.set_text(len(occurence_des_BD[list(occurence_des_BD.keys())[0]]))
#
# # Affichez le diagramme de Venn
# plt.show()


