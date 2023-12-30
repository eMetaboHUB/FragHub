from msp_parsers import *
from set_parameters import *
from fraghubid_generator import *
import numpy as np
import time
import os

from diff_fraghubid_originals_vs_clean import *

ordered_columns = ["FILENAME",
                   "PREDICTED",
                   "FRAGHUBID",
                   "SPECTRUMID",
                   "RESOLUTION",
                   "SYNON",
                   "CHARGE",
                   "IONIZATION",
                   "MSLEVEL",
                   "FRAGMENTATIONMODE",
                   "NAME",
                   "PRECURSORMZ",
                   "EXACTMASS",
                   "AVERAGEMASS",
                   "PRECURSORTYPE",
                   "INSTRUMENTTYPE",
                   "INSTRUMENT",
                   "SMILES",
                   "INCHI",
                   "INCHIKEY",
                   "COLLISIONENERGY",
                   "FORMULA",
                   "RETENTIONTIME",
                   "IONMODE",
                   "COMMENT",
                   "NUM PEAKS",
                   "peak_list"]

# Execution de la fonction
build_window()

start_time = time.time()

# generate_fraghub_id(r"../INPUT/MSP") # one time !!!

# generating filepath list
spectrum_path_list = []
MSP_path = os.path.abspath(r"../INPUT/MSP")
for files in os.listdir(MSP_path):
    if files.endswith(".msp"):
        spectrum_path_list.append(os.path.join(MSP_path,files))

spectrum_path_list = [spectrums for spectrums in spectrum_path_list if spectrums != None]

for path in spectrum_path_list:
    spectrum_list = load_spectrum_list(path)
    break

spectrum_list = msp_cleaning_processing(spectrum_list)

# print(len(spectrum_list))

df = pd.DataFrame(spectrum_list)

df = df[ordered_columns]

print("-- MOLS HARMONIZATION AND MASS CALCULATION --")
time.sleep(0.01)
df = mols_derivation_and_calculation(df)

# print("-- MASS CALCULATION --")
# time.sleep(0.01)
# df = mass_calculator(df)

fraghub_id_column = pd.DataFrame(columns=['FRAGHUBID'])
fraghub_id_column = df[["FRAGHUBID"]]

# df.to_excel(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_msms_pos.xlsx",index=False)

print("writting soon")

chunk_size = 5000  # Taille de chaque fraction
num_chunks = int(np.ceil(df.shape[0] / chunk_size))  # Calculer le nombre de fractions

with tqdm(total=num_chunks, unit=" row", colour="green", desc="\t    writting") as pbar:
    for start in range(0, df.shape[0], chunk_size):
        df_slice = df[start:start + chunk_size]
        if start == 0:
            # Écrire les en-têtes pour la première fraction
            df_slice.to_csv(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_gnps.csv", mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)
        else:
            # Append dans le fichier sans écrire les en-têtes pour les autres fractions
            df_slice.to_csv(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_gnps.csv", mode='a', header=False, index=False, sep=";", quotechar='"', encoding="UTF-8")

        # Mettre à jour la barre de progression
        pbar.update()

# compteur = 1
# for spectrum in spectrum_list:
#     print(spectrum)
#     spectrum.to_excel(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_{compteur}.xlsx",index=False)
#     compteur += 1
#     break

diff_fraghubid_originals_vs_clean(fraghub_id_column)

print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))