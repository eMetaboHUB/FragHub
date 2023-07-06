from tqdm import tqdm
import csv
import os
import re

def complete_names(filename):
    # Dictionnaire pour stocker les InChIKey et les noms correspondants
    inchikey_names = {}

    # Lecture du fichier MSP
    with open(filename, 'r', encoding="UTF-8") as file:
        data = file.read()

    # Séparation des spectres
    spetcrum_list = data.split("\n\n")

    for spectrum in tqdm(spetcrum_list, total=len(spetcrum_list), unit=" spectrums", colour="green", desc="\t processing"):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        if re.search("\nNAME: (.*)\n",spectrum):
            name = re.search("\nNAME: (.*)\n",spectrum).group(1)

        # Si l'InChIKey et le nom sont valides, les ajouter au dictionnaire
        if inchikey and name != "None":
            inchikey_names[inchikey] = name

    # Mise à jour des noms manquants avec les noms correspondants dans la liste de dictionnaires
    updated_spetcrum_list = []
    for spectrum in tqdm(spetcrum_list, total=len(spetcrum_list), unit=" spectrums", colour="green", desc="\t updating"):
        if re.search("INCHIKEY: (.*)\n",spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n",spectrum).group(1)
        if re.search("\nNAME: (.*)\n",spectrum):
            name = re.search("\nNAME: (.*)\n",spectrum).group(1)
        if name == "None":
            if inchikey in inchikey_names.keys():
                spectrum = re.sub("\nNAME: (.*)\n",f"\nNAME: {inchikey_names[inchikey]}\n",spectrum)

        updated_spetcrum_list.append(spectrum)


    # Écriture des spectres mis à jour dans un nouveau fichier MSP
    with open(os.path.join("../names_completion/out", os.path.basename(filename).replace(".msp", "_2.msp")), 'w', encoding="UTF-8") as file:
        file.write("\n\n".join(updated_spetcrum_list))

    # Écritdes spectres mis à jour dans un fichier CSV
    with open(os.path.join("../names_completion/out", os.path.basename(filename).replace(".msp", "_2.csv")), 'w', newline='', encoding="UTF-8") as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(["SPECTRUM"])
        writer.writerows([[spectrum] for spectrum in updated_spetcrum_list])

    print("Complétion des noms terminée.")

# Exemple d'utilisation
folder = "../names_completion"

for file in os.listdir(folder):
    if file.endswith(".msp"):
        file_path = os.path.join(folder, file)

        complete_names(os.path.join(folder, file))
